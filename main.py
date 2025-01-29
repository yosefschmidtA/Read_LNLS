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
        intensidade_calculada = subset['intensitycal'].values
        intensidade_experimental = subset['intensityexp'].values

        # Acumular todas as intensidades para o cálculo do Rtotal
        todas_intensidades_exp.extend(intensidade_experimental)
        todas_intensidades_calc.extend(intensidade_calculada)

        # Normalização por theta
        somatorio1 = np.sum(np.abs(intensidade_experimental))
        somatorio2 = np.sum(np.abs(intensidade_calculada))

        if somatorio1 == 0 or somatorio2 == 0:
            continue

        nova_coluna_exp = intensidade_experimental / somatorio1
        nova_coluna_calc = intensidade_calculada / somatorio2

        resultado_final1 = np.sum((nova_coluna_exp - nova_coluna_calc) ** 2)
        resultado_final2 = np.sum(nova_coluna_exp ** 2 + nova_coluna_calc ** 2)

        resultado_final3 = resultado_final1 / resultado_final2
        resultado_final4 += resultado_final3

        resultados.append(resultado_final3)
        angulos.append(theta)

        n_total += 1

    # Cálculo do R-factor médio
    r_factor_medio = resultado_final4 / n_total if n_total > 0 else None

    # Cálculo do Rtotal usando todas as intensidades acumuladas
    todas_intensidades_exp = np.array(todas_intensidades_exp)
    todas_intensidades_calc = np.array(todas_intensidades_calc)

    somatorio1_total = np.sum(np.abs(todas_intensidades_exp))
    somatorio2_total = np.sum(np.abs(todas_intensidades_calc))

    if somatorio1_total > 0 and somatorio2_total > 0:
        nova_coluna_exp_total = todas_intensidades_exp / somatorio1_total
        nova_coluna_calc_total = todas_intensidades_calc / somatorio2_total

        resultado_final1_total = np.sum((nova_coluna_exp_total - nova_coluna_calc_total) ** 2)
        resultado_final2_total = np.sum(nova_coluna_exp_total ** 2 + nova_coluna_calc_total ** 2)

        r_factor_total = resultado_final1_total / resultado_final2_total
    else:
        r_factor_total = None

    return resultados, angulos, r_factor_medio, r_factor_total

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None

    for i in range(26, len(lines)):
        line = lines[i].strip()

        # Parar a leitura ao encontrar "fitted parameters ("
        if "fitted parameters (" in line:
            break

        if line:
            parts = line.split()

            if len(parts) == 7:
                theta_value = float(parts[3])

            elif len(parts) == 5 and theta_value is not None:
                phi = float(parts[0])
                intensitycal = float(parts[3])
                intensityexp = float(parts[4])
                data.append([phi, intensitycal, theta_value, intensityexp, True])  # Marcar como original

    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)

    # Verificar o intervalo de Phi
    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    phi_interval = phi_max - phi_min

    if phi_interval < 360 and df['Phi'].max() < 360:
        df_360 = df[df['Phi'] == 0].copy()
        df_360['Phi'] = 360
        df = pd.concat([df, df_360], ignore_index=True)

    if phi_interval == 120:
        # Replicação dos dados para cobrir 360 graus
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False, as_index=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original
        # Adicionar os novos valores ao DataFrame
        df = pd.concat([df, last_values], ignore_index=True)

        df_0_120 = df.copy()
        df_0_120['Phi'] = 120 + df_0_120['Phi']
        df_0_120['isOriginal'] = False

        df_240_360 = df.copy()
        df_240_360['Phi'] = 240 + df_240_360['Phi']

        df = pd.concat([df, df_0_120, df_240_360]).reset_index(drop=True)

    if phi_interval == 90:
        # Replicação dos dados para cobrir 360 graus
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False, as_index=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original

        # Adicionar os novos valores ao DataFrame
        df = pd.concat([df, last_values], ignore_index=True)
        df_0_90 = df.copy()
        df_0_90['Phi'] = 90 + df_0_90['Phi']

        df_90_180 = df.copy()
        df_90_180['Phi'] = 180 + df_90_180['Phi']

        df_180_270 = pd.concat([df,df_0_90]).reset_index(drop=True)
        df_180_270['Phi'] = 180 + df_180_270['Phi']

        df = pd.concat([df, df_0_90, df_180_270]).reset_index(drop=True)

    return df, r_factor_total  # Retorna r_factor_total

def interpolate_data(df, resolution=1000):
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['intensitycal']

    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')

    return phi_grid, theta_grid, intensity_grid


def plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=None, save_path=None):
    # Interpolar os dados
    plt.ion()
    phi_grid, theta_grid, intensity_grid = interpolate_data(df, resolution)

    # Criando o gráfico polar
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8), dpi=150)

    # Plotando a intensidade interpolada
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')

    # Definir o limite máximo do eixo theta com base no maior valor de theta nos dados
    max_theta = df['Theta'].max()  # Maior valor de theta presente nos dados
    ax.set_ylim(0, np.radians(max_theta))  # Limitar o eixo radial até o maior valor de theta

    # Adiciona rótulos para os ângulos theta, ajustados conforme o máximo de theta nos dados
    theta_ticks = np.linspace(0, max_theta, num=6)  # Definir até 6 ticks no eixo theta
    ax.set_yticks(np.radians(theta_ticks))  # Converte para radianos
    ax.set_yticklabels([f'{int(tick)}°' for tick in theta_ticks])  # Exibe como graus

    # Adicionando a barra de cores
    cbar = fig.colorbar(c, ax=ax, label=r'$\chi = \frac{I_{\text{exp}} - I_0}{I_0}$')
    cbar.set_label(r'$\chi(θ,φ) = \frac{I(θ,φ) - I_0(θ,φ)}{I_0(θ,φ)}$', fontsize=22)

    # Adiciona a variável fora do gráfico, no canto inferior direito
    if my_variable is not None:
        # Aqui estamos usando coordenadas relativas à figura (0 a 1)
        fig.text(0.75, 0.07, f'R-factor: {my_variable}', fontsize=22, color='black', ha='right', va='bottom',
                 fontweight='bold')

    plt.xticks(fontsize=11, fontweight='bold')
    plt.yticks(fontsize=0, fontweight='bold')
    plt.draw()

    # Salvar a figura, se o caminho de salvamento for fornecido
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')  # Salva no caminho especificado

    # Exibir a figura por 600 segundos
    plt.pause(600)


# Caminho do arquivo
file_path = 'bestBCC.out'
save_path = 'grafico_polar.png'
df, r_factor_total = process_file(file_path)
my_variable = "{:.3f}".format(r_factor_total)
plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=my_variable, save_path=save_path)
