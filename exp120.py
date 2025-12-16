from scipy.ndimage import gaussian_filter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def calcular_r_factor(df):
    resultados = []
    angulos = []
    resultado_final4 = 0
    n_total = 0

    angulos_unicos = df['Theta'].unique()

    # Acumuladores para o cálculo do Rtotal
    todas_intensidades_exp = []
    todas_intensidades_calc = []

    for theta in angulos_unicos:
        subset = df[df['Theta'] == theta]
        phi = subset['Phi'].values
        intensidade_calculada = np.round(subset['intensitycal'].values, 10)
        intensidade_experimental = np.round(subset['intensityexp'].values, 10)

        # Acumular todas as intensidades para o cálculo do Rtotal
        todas_intensidades_exp.extend(intensidade_experimental)
        todas_intensidades_calc.extend(intensidade_calculada)

        # Normalização por theta
        somatorio1 = np.round(np.sum(np.abs(intensidade_experimental)),10)
        somatorio2 = np.round(np.sum(np.abs(intensidade_calculada)),10)

        if somatorio1 == 0 or somatorio2 == 0:
            continue

        nova_coluna_exp = np.round((intensidade_experimental / somatorio1),10)
        nova_coluna_calc = np.round((intensidade_calculada / somatorio2),10)

        resultado_final1 = np.round(np.sum((nova_coluna_exp - nova_coluna_calc) ** 2),10)
        resultado_final2 = np.round(np.sum(nova_coluna_exp ** 2 + nova_coluna_calc ** 2),10)

        resultado_final3 = np.round((resultado_final1 / resultado_final2),10)
        resultado_final4 += resultado_final3

        resultados.append(resultado_final3)
        angulos.append(theta)

        n_total += 1

    # Cálculo do R-factor médio
    r_factor_medio = np.round((resultado_final4 / n_total),10) if n_total > 0 else None

    # Cálculo do Rtotal usando todas as intensidades acumuladas
    todas_intensidades_exp = np.round(np.array(todas_intensidades_exp),10)
    todas_intensidades_calc = np.round(np.array(todas_intensidades_calc),10)

    somatorio1_total = np.round(np.sum(np.abs(todas_intensidades_exp)),10)
    somatorio2_total = np.round(np.sum(np.abs(todas_intensidades_calc)),10)

    if somatorio1_total > 0 and somatorio2_total > 0:
        nova_coluna_exp_total = np.round((todas_intensidades_exp / somatorio1_total),10)
        nova_coluna_calc_total = np.round((todas_intensidades_calc / somatorio2_total),10)

        resultado_final1_total = np.round(np.sum((nova_coluna_exp_total - nova_coluna_calc_total) ** 2),10)
        resultado_final2_total = np.round(np.sum(nova_coluna_exp_total ** 2 + nova_coluna_calc_total ** 2),10)

        r_factor_total = np.round((resultado_final1_total / resultado_final2_total),10)
    else:
        r_factor_total = None

    return resultados, angulos, r_factor_medio, r_factor_total


def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None
    
    # Varre as linhas procurando o padrão do arquivo Agclean.txt
    for line in lines:
        line = line.strip()
        if not line: continue
        
        parts = line.split()

        # Identifica a linha de cabeçalho do âddngulo Theta
        # Ex: 1 40 16.6600 18.0000 1.00000 0.00000 (6 colunas)
        # O Theta é a 4ª coluna (índice 3) -> 18.0000
        if len(parts) == 6:
            try:
                # Verificamos se é numérico para evitar pegar texto errado
                possible_theta = float(parts[3])
                # Confirmação simples baseada no valor (evita pegar índices 1, 2, 3...)
                if possible_theta > 0: 
                    theta_value = possible_theta
            except ValueError:
                continue

        # Identifica a linha de dados experimentais
        # Ex: 111.000 114644. 117252. -0.0222434 (4 colunas)
        # Colunas: Phi, IntRaw, Back, Chi (Anisotropia)
        elif len(parts) == 4 and theta_value is not None:
            try:
                phi = float(parts[0])
                chi_exp = float(parts[3]) # Pegamos a última coluna que é o Chi
                
                # IMPORTANTE: Jogamos o chi_exp na coluna 'intensitycal' 
                # para reaproveitar a função de plot sem alterar tudo.
                data.append([phi, chi_exp, theta_value, chi_exp, True])
            except ValueError:
                continue

    # Cria o DataFrame
    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)
    for theta, r_factor in zip(angulos, resultados):
        print(f'Theta: {theta}, R-factor: {r_factor}')

    print(f'R-factor médio: {r_factor_medio}')
    print(f'R-factor total: {r_factor_total}')
    # Verificar o intervalo de Phi
    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    phi_interval = phi_max - phi_min
    step = None 

    # 1. Identificação do Passo de Simetria
    if phi_interval > 350:
        print("Varredura Completa. Apenas fechando bordas.")
        step = 360 # O passo é o próprio círculo
    elif 170 < phi_interval < 190:
        print("Simetria C2v (180°).")
        step = 180
    elif 80 < phi_interval < 135:
        step = 120 if phi_interval > 105 else 90
        print(f"Simetria Setorial ({step}°).")

    # 2. Replicação Robusta (Cobre buracos grandes)
    if step is not None and step < 360:
        dfs_to_concat = [df]
        current_shift = step
        # Replica até passar de 360
        while current_shift <= 360 + step:
            df_shift = df.copy()
            df_shift['Phi'] = df_shift['Phi'] + current_shift
            dfs_to_concat.append(df_shift)
            current_shift += step
        
        # Replica para trás também (importante para o 0 ser contínuo)
        df_minus = df.copy()
        df_minus['Phi'] = df_minus['Phi'] - step
        dfs_to_concat.append(df_minus)
        
        df = pd.concat(dfs_to_concat, ignore_index=True)

    # 3. "GAMBIARRA" AUTOMATIZADA (GHOST POINTS)
    # Aqui resolvemos a linha. Criamos uma margem de segurança.
    # Pegamos tudo que está perto de 360 e jogamos para o 0 (negativo)
    # Pegamos tudo que está perto de 0 e jogamos para o 360 (positivo)
    
    margin = 40 # Graus de margem de segurança
    
    # Pega o final do círculo e joga para trás do zero (ex: 350 vira -10)
    subset_end = df[df['Phi'] > (360 - margin)].copy()
    subset_end['Phi'] = subset_end['Phi'] - 360
    
    # Pega o início do círculo e joga para frente do 360 (ex: 10 vira 370)
    subset_start = df[df['Phi'] < margin].copy()
    subset_start['Phi'] = subset_start['Phi'] + 360
    
    # Junta tudo. Agora temos dados de -40 até 400 graus!
    df = pd.concat([df, subset_end, subset_start], ignore_index=True)

    # O interpolador (griddata) agora vai ter vizinhos de sobra para calcular
    # o 0 e o 360 sem deixar cicatriz nenhuma.
    return df, r_factor_total  # Retorna r_factor_total

def interpolate_data(df, resolution=1000):
    # Converte os dados de entrada para radianos
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['intensitycal']

    # --- A MUDANÇA ESTÁ AQUI ---
    # Em vez de usar min(phi) e max(phi), forçamos o grid a ser estritamente 0 a 2pi.
    # endpoint=False é importante para não duplicar o 360 agora (faremos a costura depois)
    phi_grid = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    
    # Theta continua pegando do mínimo ao máximo dos dados
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    # O griddata vai olhar para seus dados extendidos (-40 a 400) 
    # e preencher esse grid (0 a 360) perfeitamente.
    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')

    return phi_grid, theta_grid, intensity_grid
def plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=None, save_path=None):
    plt.ion()
    
    # 1. Interpolação (Agora gera apenas 0 a 360, limpo)
    phi_grid, theta_grid, intensity_grid = interpolate_data(df, resolution)
    
    sigma = 4 # Pode ajustar esse valor

    # 2. Tratamento de NaNs antes do filtro (Segurança extra)
    # Se houver algum buraco, preenchemos com 0 para o filtro não estourar
    intensity_grid = np.nan_to_num(intensity_grid, nan=0.0)

    # 3. Suavização
    if sigma > 0:
        # mode='wrap' garante que o filtro entenda que o lado esquerdo conecta com o direito
        intensity_grid = gaussian_filter(intensity_grid, sigma=sigma, mode='wrap')

    # 4. Costura Final (Stitching) para fechar o anel visualmente
    # Adicionamos o ponto 360 igual ao ponto 0
    phi_col0 = phi_grid[:, 0:1] + 2 * np.pi # Pega o início e projeta em 360
    theta_col0 = theta_grid[:, 0:1]
    inten_col0 = intensity_grid[:, 0:1]
    
    phi_grid = np.hstack([phi_grid, phi_col0])
    theta_grid = np.hstack([theta_grid, theta_col0])
    intensity_grid = np.hstack([intensity_grid, inten_col0])

    # Plotagem
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8), dpi=100)
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')

    # ... (O RESTO DO CÓDIGO PERMANECE IGUAL) ...
    max_theta = df['Theta'].max()
    ax.set_ylim(0, np.radians(max_theta))

    theta_ticks = np.linspace(0, max_theta, num=6)
    ax.set_yticks(np.radians(theta_ticks))
    ax.set_yticklabels([f'{int(tick)}°' for tick in theta_ticks])

    cbar = fig.colorbar(c, ax=ax, label='', pad=0.08)
    cbar.set_label('', fontsize=22, fontweight='bold')
    cbar.ax.yaxis.set_tick_params(labelsize=30)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    if my_variable is not None:
        fig.text(0.87, 0.03, f'{my_variable}', fontsize=34, color='black', ha='right', va='bottom', fontweight='bold')
    
    fig.text(0.94, 0.9, "Anisotropy", fontsize=34, color='black', ha='right', va='bottom', fontweight='bold')

    phi_ticks = np.linspace(0, 2 * np.pi, num=9)[:-1]
    phi_labels = [f'{int(np.degrees(tick))}°' for tick in phi_ticks]
    ax.set_xticks(phi_ticks)
    ax.set_xticklabels(phi_labels, fontsize=26, fontweight='bold')

    pad_values = [1, -1, 3, 0, -7, -6, -1, -6]
    for label, pad in zip(ax.get_xticklabels(), pad_values):
        label.set_y(label.get_position()[1] + pad * 0.01)

    ax.tick_params(pad=8)
    plt.yticks(fontsize=0, fontweight='bold')
    
    plt.draw()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')

    plt.pause(600)
# Caminho do arquivo
file_path = 'teste.txt'
save_path = 'grafico_polar3.png'
df,r_total = process_file(file_path)
#r_factor_total=0.276
my_variable = "Experimental"
plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=my_variable, save_path=save_path)
