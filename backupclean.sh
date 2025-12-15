#!/bin/bash

# 1. Verifica argumentos
if [ "$#" -ne 2 ]; then
    echo "Uso correto: ./plotR.sh <ANGULO> <ARQUIVO>"
    echo "Exemplo: ./plotR.sh 120 saida1.txt"
    exit 1
fi

ANGULO=$1
ARQUIVO=$2
SCRIPT_NAME="temp_script.py"

echo "=========================================="
echo "Configuração detectada:"
echo "Ângulo: $ANGULO"
echo "Arquivo: $ARQUIVO"
echo "=========================================="

# 2. Decide qual script gerar
if [ "$ANGULO" -ge 115 ] && [ "$ANGULO" -le 125 ]; then
    echo "-> Selecionado modo TRIANGULAR (aprox. 120 graus)."
    echo "-> Gerando script baseado em new.py..."
    
    # ATENÇÃO: O código abaixo usa espaços estritos (sem tabs)
    cat << 'EOF' > $SCRIPT_NAME
import sys
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
    todas_intensidades_exp = []
    todas_intensidades_calc = []

    for theta in angulos_unicos:
        subset = df[df['Theta'] == theta]
        phi = subset['Phi'].values
        intensidade_calculada = np.round(subset['intensitycal'].values, 10)
        intensidade_experimental = np.round(subset['intensityexp'].values, 10)

        todas_intensidades_exp.extend(intensidade_experimental)
        todas_intensidades_calc.extend(intensidade_calculada)

        somatorio1 = np.round(np.sum(np.abs(intensidade_experimental)), 10)
        somatorio2 = np.round(np.sum(np.abs(intensidade_calculada)), 10)

        if somatorio1 == 0 or somatorio2 == 0:
            continue

        nova_coluna_exp = np.round((intensidade_experimental / somatorio1), 10)
        nova_coluna_calc = np.round((intensidade_calculada / somatorio2), 10)

        resultado_final1 = np.round(np.sum((nova_coluna_exp - nova_coluna_calc) ** 2), 10)
        resultado_final2 = np.round(np.sum(nova_coluna_exp ** 2 + nova_coluna_calc ** 2), 10)

        resultado_final3 = np.round((resultado_final1 / resultado_final2), 10)
        resultado_final4 += resultado_final3

        resultados.append(resultado_final3)
        angulos.append(theta)
        n_total += 1

    r_factor_medio = np.round((resultado_final4 / n_total), 10) if n_total > 0 else None

    todas_intensidades_exp = np.round(np.array(todas_intensidades_exp), 10)
    todas_intensidades_calc = np.round(np.array(todas_intensidades_calc), 10)

    somatorio1_total = np.round(np.sum(np.abs(todas_intensidades_exp)), 10)
    somatorio2_total = np.round(np.sum(np.abs(todas_intensidades_calc)), 10)

    if somatorio1_total > 0 and somatorio2_total > 0:
        nova_coluna_exp_total = np.round((todas_intensidades_exp / somatorio1_total), 10)
        nova_coluna_calc_total = np.round((todas_intensidades_calc / somatorio2_total), 10)
        
        r1_total = np.round(np.sum((nova_coluna_exp_total - nova_coluna_calc_total) ** 2), 10)
        r2_total = np.round(np.sum(nova_coluna_exp_total ** 2 + nova_coluna_calc_total ** 2), 10)
        r_factor_total = np.round((r1_total / r2_total), 10)
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
        if "fitted parameters (" in line:
            break

        if line:
            parts = line.split()
            if len(parts) == 7:
                theta_value = float(parts[3])
            elif len(parts) == 5 and theta_value is not None:
                phi = float(parts[0])
                intensi1 = float(parts[1])
                intensi2 = float(parts[2])
                # Evita divisao por zero se intensi2 for 0
                intensitycal = float((intensi1 - intensi2) / intensi2) if intensi2 != 0 else 0
                intensityexp = float(parts[4])
                data.append([phi, intensitycal, theta_value, intensityexp, True])

    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)
    
    for theta, r_factor in zip(angulos, resultados):
        print(f'Theta: {theta}, R-factor: {r_factor}')

    print(f'R-factor médio: {r_factor_medio}')
    print(f'R-factor total: {r_factor_total}')

    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    phi_interval = phi_max - phi_min

    # Lógica de Replicação e Ajuste de Intervalo
    # Se intervalo for pequeno (ex: 117 ou 120), replicamos
    if 110 < phi_interval < 125:
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False, as_index=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()
        last_values['Phi'] = first_values['Phi']
        df = pd.concat([df, last_values], ignore_index=True)

        df_0_120 = df.copy()
        df_0_120['Phi'] = 120 + df_0_120['Phi']
        df_0_120['isOriginal'] = False

        df_240_360 = df.copy()
        # Ajuste fino para o caso 117 vs 120
        offset_third = 243 if phi_interval == 117 else 240
        df_240_360['Phi'] = offset_third + df_240_360['Phi']

        df = pd.concat([df, df_0_120, df_240_360]).reset_index(drop=True)

    # Adicione aqui as outras lógicas (90 graus) se necessário
    
    return df, r_factor_total

def interpolate_data(df, resolution=1000):
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['intensitycal']

    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)
    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)
    
    # Aqui usamos cubic conforme seu ultimo pedido
    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')
    return phi_grid, theta_grid, intensity_grid

def plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=None, save_path=None):
    plt.ion()
    # Verifica e aplica simetria se necessário (código simplificado para o caso 120 já replicado)
    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    intervalo = phi_max - phi_min
    
    if intervalo > 300:
        df_plot = df.copy()
    else:
        # Fallback para espelhamento simples se não tiver replicado
        df_left = df[(df['Phi'] >= 90) & (df['Phi'] <= 270)].copy()
        df_right = df_left.copy()
        df_right['Phi'] = (df_right['Phi'] + 180) % 360
        df_plot = pd.concat([df_left, df_right], ignore_index=True)

    # Tapa-buraco (0 e 360)
    if not df_plot.empty:
        min_phi = df_plot['Phi'].min()
        subset_start = df_plot[df_plot['Phi'] == min_phi].copy()
        subset_start['Phi'] = subset_start['Phi'] + 360
        max_phi = df_plot['Phi'].max()
        subset_end = df_plot[df_plot['Phi'] == max_phi].copy()
        subset_end['Phi'] = subset_end['Phi'] - 360
        df_plot = pd.concat([df_plot, subset_start, subset_end], ignore_index=True)

    phi_grid, theta_grid, intensity_grid = interpolate_data(df_plot, resolution)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8), dpi=100)
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')

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
        fig.text(0.87, 0.03, f'R-factor: {my_variable}', fontsize=34, color='black', ha='right', va='bottom', fontweight='bold')
    
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
    
    # Pausa reduzida para teste
    plt.pause(10)

# EXECUÇÃO PRINCIPAL
file_path = sys.argv[1]
save_path = 'grafico_polar_120.png'
df, r_factor_total = process_file(file_path)
my_variable = "{:.3f}".format(r_factor_total) if r_factor_total is not None else "N/A"
plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=my_variable, save_path=save_path)
EOF

else
    echo "-> Selecionado modo PADRÃO (360, 90, 180 graus)."
    echo "-> Gerando script padrão..."
    
    # ATENÇÃO: Você precisa colocar o código do plotR.py aqui!
    # Vou colocar um placeholder para não dar erro de sintaxe, 
    # mas VOCÊ DEVE SUBSTITUIR PELO SEU CÓDIGO REAL.
    cat << 'EOF' > $SCRIPT_NAME
import sys
print("ERRO: O código para o modo padrão (else) ainda não foi colado no plotR.sh!")
print("Por favor, edite o plotR.sh e cole o conteúdo de plotR.py no bloco else.")
EOF

fi

# 3. Executa
echo "-> Executando o Python..."
python3 $SCRIPT_NAME "$ARQUIVO"

# Limpa
rm $SCRIPT_NAME
echo "-> Concluído!"
