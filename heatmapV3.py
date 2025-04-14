import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import os
import matplotlib.ticker as ticker

interlayer1 = []
interlayer2 = []
rfactor = []

with open('dados.txt', 'r') as file:
    lines = file.readlines()
    for line in lines:
        # Divide cada linha em partes
        partes = line.split()
        if len(partes) == 3:
            # Lê os valores de interlayer1, interlayer2 e rfactor
            interlayer1_value = float(partes[0])
            interlayer2_value = float(partes[1])
            rfactor_value = float(partes[2])

            # Adiciona os valores às listas
            interlayer1.append(interlayer1_value)
            interlayer2.append(interlayer2_value)
            rfactor.append(rfactor_value)

# Converte as listas em arrays numpy para facilitar o manuseio
interlayer1 = np.array(interlayer1)  # Multiplicando por 3.905 conforme no seu código original
interlayer2 = np.array(interlayer2)
rfactor = np.array(rfactor)

# Cria uma grade para interpolar os dados
xi = np.linspace(min(interlayer1), max(interlayer1), 100)
yi = np.linspace(min(interlayer2), max(interlayer2), 100)
xi, yi = np.meshgrid(xi, yi)
zi = griddata((interlayer1, interlayer2), rfactor, (xi, yi), method='linear')

# Define os níveis dos contornos
num_levels = 100  # Reduz o número de níveis de contorno
levels = np.linspace(min(rfactor), max(rfactor), num_levels)

# Cria o gráfico de contorno preenchido
plt.figure(figsize=(10, 5))
contour_filled = plt.contourf(xi, yi, zi, levels=levels, cmap='jet_r')

# Adiciona as linhas de contorno
contour_lines = plt.contour(xi, yi, zi, levels=levels, colors='white', linewidths=0.0)

# Adiciona a barra de cores
cbar = plt.colorbar(contour_filled)
cbar.set_label('R$_{factor}$', rotation=0, labelpad=0, fontsize=30)
cbar.ax.yaxis.set_label_position('left')  # Move the colorbar label to the left
cbar.ax.yaxis.set_label_coords(0.4, 1.01)  # Ajusta a posição do rótulo acima da barra de cores
cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
cbar.ax.tick_params(labelsize=28)

#Vitor, pra fazer os gráficos para outros layers além do d12d23, lembra de mudar aqui os nomes
plt.xlabel('D$_{1-2}$ (Å)', fontsize=46)
plt.ylabel('D$_{2-3}$ (Å)', fontsize=46)
save_dir = "relaxation"
os.makedirs(save_dir, exist_ok=True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
filename = f"relax0.jpg"
filepath = os.path.join(save_dir, filename)
# Define o mínimo global do R-factor
rfactor_min = np.min(rfactor)


# Também marca o ponto exato do mínimo com um ponto preto
min_index = np.argmin(rfactor)
plt.plot(interlayer1[min_index], interlayer2[min_index], 'ko', markersize=8, label='Global minimum')

# Adiciona legenda para o ponto mínimo
plt.legend(fontsize=14, bbox_to_anchor=(0.95, 1.15))
plt.savefig(filepath, dpi=250, bbox_inches='tight')

# Exibe o gráfico
plt.show()
