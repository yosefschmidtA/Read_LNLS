from scipy.interpolate import griddata
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('TkAgg')
plt.ion()

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None

    for i in range(17, len(lines)):
        line = lines[i].strip()

        if line:
            parts = line.split()

            if len(parts) == 6:
                theta_value = float(parts[3])

            elif len(parts) == 4 and theta_value is not None:
                phi = float(parts[0])
                col1 = float(parts[1])
                col2 = float(parts[2])
                intensity = float(parts[3])
                data.append([phi, col1, col2, theta_value, intensity, True])  # Marcar como original

    df = pd.DataFrame(data, columns=['Phi', 'Col1', 'Col2', 'Theta', 'Intensity', 'IsOriginal'])
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
        df = df.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
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
        df = df.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
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

    return df

# Função para interpolar os dados
def interpolate_data(df, resolution=1000):
    # Definir uma grade regular para a interpolação
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['Intensity']

    # Criando uma grade de pontos onde queremos interpolar
    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    # Criando uma grade de malha para interpolação
    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    # Realizar a interpolação
    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='linear')

    return phi_grid, theta_grid, intensity_grid


# Função para gerar o gráfico polar
def plot_polar_interpolated(df, resolution=500):
    # Interpolar os dados
    plt.ion()
    phi_grid, theta_grid, intensity_grid = interpolate_data(df, resolution)
    # Criando o gráfico polar
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10,8), dpi=100)

    # Plotando a intensidade interpolada
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')


    # Definir o limite máximo do eixo theta com base no maior valor de theta nos dados
    max_theta = df['Theta'].max()  # Maior valor de theta presente nos dados
    ax.set_ylim(0, np.radians(max_theta))  # Limitar o eixo radial até o maior valor de theta

    # Adiciona rótulos para os ângulos theta, ajustados conforme o máximo de theta nos dados
    theta_ticks = np.linspace(0, max_theta, num=6)  # Definir até 6 ticks no eixo theta
    ax.set_yticks(np.radians(theta_ticks))  # Converte para radianos
    ax.set_yticklabels([f'{int(tick)}°' for tick in theta_ticks],)  # Exibe como graus


    ax.set_xlabel('θ', fontsize=36)
    ax.set_ylabel('φ', fontsize=36, labelpad=1, rotation = 360)
    ax.yaxis.set_label_coords(0.8, 0.93)  # Ajustar a posição do rótulo 'φ'

    # Adicionando a barra de cores
    cbar = fig.colorbar(c, ax=ax, label=r'$\chi = \frac{I_{\text{exp}} - I_0}{I_0}$')
    cbar.set_label(r'$\chi = \frac{I_{\text{exp}} - I_0}{I_0}$', fontsize=36)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=16)
    plt.draw()
    plt.pause(600)
    plt.close()

df = process_file("exp_Fe2P_fitted.out")
plot_polar_interpolated(df)
