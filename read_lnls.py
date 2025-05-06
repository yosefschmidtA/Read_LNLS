def carregar_config(arquivo):
    parametros = {}
    with open(arquivo, 'r') as f:
        for linha in f:
            # Remove comentários (tudo após #)
            linha = linha.split('#', 1)[0].strip()
            if '=' in linha and linha:  # Garante que a linha não está vazia
                chave, valor = linha.split('=', 1)
                chave = chave.strip()
                valor = valor.strip()
                try:
                    # Tenta converter para número (int ou float)
                    if '.' in valor:
                        parametros[chave] = float(valor)
                    else:
                        parametros[chave] = int(valor)
                except ValueError:
                    # Se não for número, mantém como string
                    parametros[chave] = valor
    return parametros

# Ler parâmetros do arquivo config.txt
parametros = carregar_config('config.txt')
import matplotlib.colors as mcolors
# Acessando os parâmetros no código
file_prefix = parametros["file_prefix"]
thetai = parametros["thetai"]
thetaf = parametros["thetaf"]
dtheta = parametros["dtheta"]
phii = parametros["phii"]
phif = parametros["phif"]
dphi = parametros["dphi"]
channel = parametros["channel"]
symmetry = parametros["symmetry"]
indice_de_plotagem = parametros["indice_de_plotagem"]
shirley_tempo = parametros["shirley_tempo"]
poli_tempo = parametros["poli_tempo"]
fft_tempo = parametros["fft_tempo"]
arquivo_saida = parametros["arquivo_saida"]
rotate_angle = parametros["rotate_angle"]

def fourier_symmetrization(theta_values, phi_values, intensity_values, symmetry):
    """
    Aplica a simetrização por expansão em Fourier nos dados XPD.

    Parameters:
        theta_values (array): Valores de theta.
        phi_values (array): Valores de phi.
        intensity_values (2D array): Intensidades organizadas como [theta][phi].
        symmetry (int): Grau de simetria (ex.: 4 para C4).

    Returns:
        intensity_symmetric (2D array): Intensidades simetrizadas.
    """
    n_theta = len(theta_values)
    n_phi = len(phi_values)

    # Inicializar array para intensidades simetrizadas
    intensity_symmetric = np.zeros_like(intensity_values)

    for i, theta in enumerate(theta_values):
        # Extrair a curva de intensidade para o theta atual
        f = intensity_values[i, :]

        # Calcular a Transformada de Fourier
        F = np.fft.fft(f)

        # Criar uma cópia para armazenar apenas os componentes simétricos
        F_symmetric = np.zeros_like(F, dtype=complex)

        # Manter apenas os harmônicos que atendem à simetria (ex.: múltiplos de symmetry)
        for u in range(0, n_phi):
            if u % symmetry == 0:
                F_symmetric[u] = F[u]

        # Calcular a Transformada Inversa de Fourier com os componentes simétricos
        f_symmetric = np.fft.ifft(F_symmetric).real

        # Salvar a curva simetrizada
        intensity_symmetric[i, :] = f_symmetric
        save_dir = "FFT"
        #os.makedirs(save_dir, exist_ok=True)
        '''# Plotar as curvas original e simetrizada (opcional)
        plt.figure()
        plt.plot(phi_values, f, label='Experimental', linestyle='--')
        plt.plot(phi_values, f_symmetric, label=f'FFT (C{symmetry})')
        plt.title(f'Theta = {theta:.2f}')
        plt.xlabel('Phi (°)')
        plt.ylabel('Intensity')
        plt.legend()

        filename = f"FFT{theta}.jpg"
        filepath = os.path.join(save_dir, filename)
        # Salvar a figura
        plt.savefig(filepath, dpi=100, bbox_inches='tight')
        plt.close()'''




    return intensity_symmetric

def gaussian_fit(x, amplitude, mean, stddev):
    """
    Retorna os valores de uma gaussiana para os parâmetros dados.
    amplitude: Altura máxima do pico.
    mean: Posição do pico.
    stddev: Largura do pico (desvio padrão).
    """
    return amplitude * np.exp(-0.5 * ((x - mean) / stddev) ** 2)
def generate_file_names(prefix, thetai, thetaf, dtheta, phii, phif, dphi):
    theta_values = [thetai + i * dtheta for i in range((thetaf - thetai) // dtheta + 1)]
    phi_values = [phii + j * dphi for j in range((phif - phii) // dphi + 1)]
    file_names = []
    for theta in theta_values:
        for phi in phi_values:
            file_name = f"{prefix}{theta}.{phi}"
            file_names.append(file_name)
    return file_names
def shirley_background(x_data, y_data, init_back, end_back, n_iterations=6):
    """
    Calcula o fundo de Shirley para um espectro de intensidade.

    Parameters:
    x_data (array): O vetor de energias de ligação (ou qualquer outra variável x).
    y_data (array): O vetor de intensidades do espectro.
    init_back (int): Índice inicial para calcular o fundo.
    end_back (int): Índice final para calcular o fundo.
    n_iterations (int): Número de iterações para refinar o fundo.

    Returns:
    background (array): O vetor de fundo calculado.
    """
    # Inicializa o vetor de fundo com os valores nas extremidades
    background = np.zeros_like(y_data)
    background0 = np.zeros_like(y_data)

    # Definindo os valores de intensidade nas extremidades
    a = y_data[init_back]
    b = y_data[end_back]

    # Calcula o fundo de Shirley por n iterações
    for nint in range(n_iterations):
        for k2 in range(end_back, init_back - 1, -1):
            sum1 = sum(y_data[k] - background0[k] for k in range(end_back, k2 - 1, -1))
            sum2 = sum(y_data[k] - background0[k] for k in range(end_back, init_back - 1, -1))

            # Calcula o fundo interpolado entre as extremidades
            if sum2 != 0:
                background[k2] = (a - b) * sum1 / sum2 + b

        # Ajuste apenas para garantir suavidade nas extremidades
        background[:init_back] = background[init_back]

        # Atualiza o fundo de referência para a próxima iteração
        background0 = background.copy()

    return background
def smooth(data, sigma=2):
    return gaussian_filter1d(data, sigma=sigma)
# Função para o polinômio de grau 3
def polynomial_3(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d


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
output_file_path = '../simetrizados.txt'  # não mexa nesse arquivo
# Função para gerar os nomes dos arquivos esperados

# Gerar os nomes de arquivos esperados
expected_files = generate_file_names(file_prefix, thetai, thetaf, dtheta, phii, phif, dphi)

# Verificar quais arquivos existem no diretório atual
existing_files = [f for f in expected_files if os.path.isfile(f)]

# Listas para armazenar os dados separados
data_one_column = []
output_file_xps = "../saidaxps.txt"

with open(output_file_xps, 'w') as log_file:
    for file in existing_files:

        try:
            # Assumindo que o formato do arquivo é 'prefixo_theta.phi'
            theta, phi = file[len(file_prefix):].split('.')
            theta = int(theta)
            phi = int(phi)
        except ValueError:
            continue  # Caso o nome do arquivo não tenha o formato esperado

        with open(file, 'r') as f:
            data_valid = False  # Flag para verificar se estamos no conjunto com o canal escolhido
            contador_banda = 0
            data_one_column = []

            for line in f:
                # Remove espaços e divide por colunas
                columns = line.strip().split()

                if len(columns) == 7:  # Linha com 7 colunas (cabeçalho)
                    if data_valid:
                        break  # Encerra a coleta ao encontrar um novo cabeçalho

                    # Verifique se a primeira coluna tem o valor do canal desejado
                    if float(columns[0]) == channel:
                        data_valid = True  # Ativa o flag para processar os dados subsequentes
                        log_file.write(f"\n{theta} {phi} {channel}\n")

                elif data_valid and len(columns) == 1:  # Linha com 1 coluna
                    value = float(columns[0])
                    contador_banda += 1
                    data_one_column.append(value)
                    log_file.write(f"{value} {contador_banda}\n")


def process_file(file_name, output_file):
    """
    Lê o arquivo de entrada e processa os dados aplicando o fundo Shirley.
    """

    def process_block(block, theta_values, phi_values):
        def doniach_sunjic(x, amp, mean, gamma, beta):
            denom = (x - mean) ** 2 + gamma ** 2
            return (amp / np.pi) * (gamma / denom) * (1 + beta * (x - mean) / denom)

        def double_doniach(x, amp1, mean1, gamma1, beta1, amp2, mean2, gamma2, beta2):
            return (doniach_sunjic(x, amp1, mean1, gamma1, beta1) + doniach_sunjic(x, amp2, mean2, gamma2, beta2))

        y_values = np.array([row[0] for row in block])  # Intensidades
        x_values = np.array([row[1] for row in block])  # Índices/canais

        y_smoothed_raw = smooth(y_values, sigma=2.0)
        init_back = 1
        end_back = len(x_values) - 1

        shirley_bg_smoothed = shirley_background(x_values, y_smoothed_raw, init_back, end_back)
        y_corrected_smoothed = y_smoothed_raw - shirley_bg_smoothed

        # Forçar valores positivos (removendo possíveis artefatos negativos)
        positive_values = y_corrected_smoothed.copy()
        positive_values[positive_values < 0] = 0

        fitted_double_doniach = fitted_double_gaussian = np.zeros_like(x_values)

        # Ajuste Doniach-Sunjic apenas se indice_de_plotagem == 1

        if indice_de_plotagem == 3:
            initial_guess_doniach = [210000, 15, 1, 0.1, 10000, 31, 1, 0.1]
            bounds_doniach = ([200000, 10, 0.5, 0, 5000, 30, 0.5, 0], [300000, 20, 8, 2, 50000, 32, 4, 2])
            try:
                popt_doniach, _ = curve_fit(double_doniach, x_values, positive_values, p0=initial_guess_doniach,
                                            bounds=bounds_doniach)
                amp1, mean1, gamma1, beta1, amp2, mean2, gamma2, beta2 = popt_doniach
                doniach1 = doniach_sunjic(x_values, amp1, mean1, gamma1, beta1)
                doniach2 = doniach_sunjic(x_values, amp2, mean2, gamma2, beta2)
                fitted_double_doniach = doniach1 + doniach2
            except Exception as e:
                print(f"Erro no ajuste das Doniach-Sunjic: {e}")

        # Ajuste Gaussianas apenas se indice_de_plotagem == 2
        if indice_de_plotagem == 2:
            def double_gaussian(x, amp1, mean1, std1, amp2, mean2, std2):
                return (gaussian_fit(x, amp1, mean1, std1) + gaussian_fit(x, amp2, mean2, std2))

            initial_guess = [5300, 9, 1,                 5500, 12.7, 1]#2  #1
            bounds = ([5000, 8, 1,       3000, 12.5, 1],
                      [100000, 10, 4,      6000, 14, 3])
                       #A(i)P(i)La(i) #A(i2)P(i2)La(i2)---
            try:
                popt, _ = curve_fit(double_gaussian, x_values, positive_values, p0=initial_guess, bounds=bounds)
                amp1, mean1, std1, amp2, mean2, std2 = popt
                gaussian1 = gaussian_fit(x_values, amp1, mean1, std1)
                gaussian2 = gaussian_fit(x_values, amp2, mean2, std2)
                fitted_double_gaussian = gaussian1 + gaussian2
            except Exception as e:
                print(f"Erro no ajuste das gaussianas: {e}")

        total_area = trapezoid(positive_values, x_values)
        print("Area: ", total_area)

        if indice_de_plotagem == 0:
            return list(zip(y_corrected_smoothed, x_values)), total_area

        # Plotagem - Mantendo a estrutura original
        if indice_de_plotagem == 3:
            title = f"XPS with Doniach-Sunjic (θ={theta_values}, φ={phi_values})"
            plt.figure(figsize=(10, 6))
            plt.plot(x_values, y_smoothed_raw, label='Original', marker='o')
            plt.plot(x_values, shirley_bg_smoothed, label='Shirley Background', linestyle='--')
            plt.plot(x_values, y_corrected_smoothed, label='Corrected', marker='x')
            plt.plot(x_values, doniach1, label='Doniach 1', linestyle='-.', color='blue')
            plt.plot(x_values, doniach2, label='Doniach 2', linestyle=':', color='green')
            plt.plot(x_values, fitted_double_doniach, label='Sum of Doniach-Sunjic functions', color='red')
            plt.fill_between(x_values, positive_values, color='yellow', alpha=0.5, label='Corrected Area')
            plt.xlabel('Channel')
            plt.ylabel('Intensity')
            plt.title(title)
            plt.legend()
            plt.grid(True)
            plt.show()
            plt.pause(shirley_tempo)
            plt.close()

        if indice_de_plotagem == 4:
            title = f"XPS with Gaussians components (θ={theta_values}, φ={phi_values})"
            plt.figure(figsize=(10, 6))
            plt.plot(x_values, y_smoothed_raw, label='Original', marker='o')
            plt.plot(x_values, shirley_bg_smoothed, label='Shirley Background', linestyle='--')
            plt.plot(x_values, y_corrected_smoothed, label='Corrected', marker='x')
            plt.plot(x_values, gaussian1, label='Gaussian 1', linestyle='-.', color='blue')
            plt.plot(x_values, gaussian2, label='Gaussian 2', linestyle=':', color='green')
            plt.plot(x_values, fitted_double_gaussian, label='Sum of the Gaussians', color='red')
            plt.fill_between(x_values, positive_values, color='yellow', alpha=0.5, label='Corrected Area')
            plt.xlabel('Channel')
            plt.ylabel('Intensity')
            plt.title(title)
            plt.legend()
            plt.grid(True)


        if indice_de_plotagem == 1:
            #plt.ion()
            save_dir = "XPS_Shirley"
            os.makedirs(save_dir, exist_ok=True)
            title = f"X-ray Photoelectron Spectrocopy: (θ={theta_values}, φ={phi_values})"
            plt.figure(figsize=(8, 6))
            plt.plot(x_values, y_smoothed_raw, label='Original', marker='o')
            plt.plot(x_values, shirley_bg_smoothed, label='Shirley Background', linestyle='--')
            plt.plot(x_values, y_corrected_smoothed, label='Correction', marker='x')
            plt.fill_between(x_values, positive_values, color='yellow', alpha=0.5, label='Corrected Area')
            plt.xlabel('Channel')
            plt.ylabel('Intensity')
            plt.title(title)
            plt.legend()
            plt.grid(True)
            #Nome do arquivo (usando os valores de θ e φ para identificação)
            filename = f"XPS_Shirley_theta_{theta_values}_phi_{phi_values}.jpg"
            filepath = os.path.join(save_dir, filename)
            # Salvar a figura
            plt.savefig(filepath, dpi=100, bbox_inches='tight')
            plt.show()


        return list(zip(y_corrected_smoothed, x_values)), total_area

    with open(file_name, 'r') as file:
        data = file.readlines()

    corrected_data = []
    block = []
    header = None

    # Processa cada linha do arquivo
    for line in data:
        columns = line.strip().split()

        if len(columns) == 3:
            # Cabeçalho do bloco
            if block:  # Se houver um bloco acumulado, processa
                corrected_data.append((header, process_block(block, theta_values, phi_values)))
                block = []
            header = line.strip()  # Salva o cabeçalho do bloco atual
        elif len(columns) == 2:
            # Adiciona valores ao bloco atual
            block.append([float(columns[0]), int(columns[1])])

        if len(columns) == 3:
            theta_values = columns[0]  # Ajuste conforme necessário
            phi_values = columns[1]   # Ajuste conforme necessário

    if block:
        corrected_data.append((header, process_block(block, theta_values, phi_values)))

    # Grava os resultados corrigidos no arquivo de saída
    with open(output_file, 'w') as out_file:
        for header, block_data in corrected_data:
            _, total_area = block_data
            out_file.write(f"{header} {total_area:.1f}\n")
            for y_corr, x in block_data[0]:
                out_file.write(f"{y_corr:.2f} {x}\n")



# Arquivos de entrada e saída
file_name = "../saidaxps.txt"
output_file = "../saidashirley.txt"

# Processa o arquivo
process_file(file_name, output_file)

#process_file_2 lê o arquivo saidashirley.txt que está com os dados após a subtração do fundo shirley e salva no dataframe "df"

def process_file_2(output_file):
    """
    Lê o arquivo de entrada, processa os dados e salva em um DataFrame.
    Modificado para salvar apenas o valor da primeira, segunda e quarta coluna.
    """
    with open(output_file, 'r') as file:
        data = file.readlines()

    # Lista para armazenar os resultados
    results = []

    # Processa cada linha do arquivo
    for line in data:
        columns = line.strip().split()

        if len(columns) == 4:  # Linha com 4 colunas (dados de theta, phi, channel e intensidade)
            # Salva o valor da primeira, segunda e quarta coluna
            theta = float(columns[0])  # Primeira coluna
            phi = float(columns[1])    # Segunda coluna
            intensity = float(columns[3])  # Quarta coluna

            # Adiciona os resultados na lista
            results.append({'theta': theta, 'phi': phi, 'intensity': intensity})

    # Cria o DataFrame com os dados
    df = pd.DataFrame(results)

    return df

# Processa o arquivo e salva os resultados em um DataFrame
df = process_file_2("../saidashirley.txt")

# Salva o DataFrame em um arquivo .txt
output_txt_file = "../saidatpintensity.txt"

with open(output_txt_file, 'w') as f:
    f.write(df.to_string(index=False))  # Salva os dados sem o índice

print(f"Dados salvos em {output_txt_file}")


def process_and_plot(input_file, output_file, plot_dir="plots", phi_values_to_evaluate=None):
    # Lê os dados do arquivo, pulando a primeira linha
    data = pd.read_csv(input_file, sep='\s+', skiprows=1, names=['theta', 'phi', 'intensity'])

    # Agrupa os dados por theta
    grouped = data.groupby('theta')

    # Lista para armazenar os resultados
    results = []

    # Processa cada grupo de theta
    for theta, group in grouped:
        phi = group['phi'].values
        intensity = group['intensity'].values

        # Ajuste polinomial
        try:
            popt, _ = curve_fit(polynomial_3, phi, intensity)
            a, b, c, d = popt
            results.append({'theta': theta, 'a': a, 'b': b, 'c': c, 'd': d})

            # Criando valores de phi de acordo com o intervalo definido
            phi_fine = np.arange(phii, phif + dphi, dphi)  # Usando phi_values
            intensity_fitted = polynomial_3(phi_fine, *popt)
            save_dir = "PoliThird"
            '''os.makedirs(save_dir, exist_ok=True)
            # Plotando os dados e o ajuste
            plt.figure(figsize=(8, 6))
            plt.plot(phi, intensity, linestyle='-', color='blue', alpha=0.5, label="Experimental")
            plt.plot(phi_fine, intensity_fitted, label="Polynomial Fitting", color="red", linewidth=2)
            plt.title(f"Third-Degree Polynomial Fit - Theta = {theta}")
            plt.xlabel("Phi")
            plt.ylabel("Intensity")
            plt.legend()
            plt.grid()'''
            #filename = f"PoliCurve{theta}.jpg"
            #filepath = os.path.join(save_dir, filename)
            # Salvar a figura
            #plt.savefig(filepath, dpi=100, bbox_inches='tight')
            #plt.show()
            #plt.pause(poli_tempo)
            #plt.close()


            # Se phi_values_to_evaluate for fornecido, calcule os valores para os pontos de phi fornecidos
            if phi_values_to_evaluate is not None:
                for phi_value in phi_values_to_evaluate:
                    intensity_at_phi = polynomial_3(phi_value, *popt)
                    print(f"Valor da intensidade para phi = {phi_value} (theta = {theta}): {intensity_at_phi}")

        except Exception as e:
            print(f"Erro ao ajustar os dados para theta = {theta}: {e}")
            continue

    # Cria um DataFrame com os coeficientes ajustados
    results_df = pd.DataFrame(results)

    # Salva os coeficientes no arquivo de saída
    results_df.to_csv(output_file, index=False, float_format='%.6f')
    print(f"Resultados salvos em {output_file}")

    # Formato do cabeçalho do arquivo de saída
    num_theta = results_df['theta'].nunique()  # Número de ângulos theta únicos
    num_phi = len(phi_fine)  # Número de ângulos phi únicos
    num_points = len(data)  # Total de pontos
    theta_initial = results_df['theta'].min()  # Valor inicial de Theta
    phi_values = np.arange(phii, phii + dphi, dphi)

    # Salvando os dados ajustados no formato esperado
    with open(output_file, 'w') as file:
        # Cabeçalho inicial
        file.write(f"      {num_theta}    {num_points}    0     datakind beginning-row linenumbers\n")
        file.write(f"----------------------------------------------------------------\n")
        file.write(f"MSCD Version 1.00 Yufeng Chen and Michel A Van Hove\n")
        file.write(f"Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720\n")
        file.write(f"Copyright (c) Van Hove Group 1997. All rights reserved\n")
        file.write(f"--------------------------------------------------------------\n")
        file.write(f"angle-resolved photoemission extended fine structure (ARPEFS)\n")
        file.write(f"experimental data for Fe 2p3/2 from Fe on STO(100)  excited with hv=1810eV\n")
        file.write(f"\n")
        file.write(f"provided by Pancotti et al. (LNLS in 9, June 2010)\n")
        file.write(f"   intial angular momentum (l) = 1\n")
        file.write(f"   photon polarization angle (polar,azimuth) = (  30.0,   0.0 ) (deg)\n")
        file.write(f"   sample temperature = 300 K\n")
        file.write(f"\n")
        file.write(f"   photoemission angular scan curves\n")
        file.write(f"     (curve point theta phi weightc weighte//k intensity chiexp)\n")
        file.write(f"      {num_theta}     {num_points}       1       {num_theta}     {num_phi}     {num_points}\n")

        # Loop para os diferentes valores de θ
        number_of_theta = 0
        for theta in sorted(results_df['theta'].unique()):
            number_of_theta += 1
            first_row = results_df[results_df['theta'] == theta].iloc[0]
            file.write(
                f"       {number_of_theta}     {num_phi}       19.5900     {first_row['theta']:.4f}      1.00000      0.00000\n"
            )
            subset = results_df[results_df['theta'] == theta]

            # Obtendo os dados experimentais de intensity e phi
            group = data[data['theta'] == theta]
            phi = group['phi'].values
            intensity = group['intensity'].values

            # Recalcula phi_fine e intensity_fitted aqui dentro do loop
            phi_fine = np.arange(phii, phif + dphi, dphi)

            for _, row in subset.iterrows():
                # Pega os coeficientes ajustados do polinômio
                a, b, c, d = row['a'], row['b'], row['c'], row['d']

                # Calcular intensidade ajustada para os valores de phi
                intensity_fitted = polynomial_3(phi_fine, a, b, c, d)

                # Calcular a média da intensidade ajustada
                mean_intensity = np.mean(intensity)
                mean_intensity2 = np.mean(intensity_fitted)
                # Calcular Chi para cada valor de phi
                Chi = ((intensity - intensity_fitted) / intensity_fitted)
                Chi2 =((intensity-mean_intensity) / mean_intensity)
                Chi3 = ((intensity-mean_intensity2) / mean_intensity2)

            # Escreve cada valor de phi_fine, intensity_fitted, mean_intensity e Chi em uma linha separada
            for p, i, chi in zip(phi, intensity, Chi):
                file.write(f"      {p:.5f}      {i:.1f}      {mean_intensity:.1f}      {chi:.7f}\n")
            file.write("")  # Linha em branco para separar os blocos de dados


# Arquivos de entrada e saída
input_file = "../saidatpintensity.txt"  # Substitua pelo nome do seu arquivo
output_file = "../coeficientes_ajustados.txt"
plot_dir = "../plots"

os.makedirs(plot_dir, exist_ok=True)

phi_values = np.arange(phii, phii + dphi, dphi)


# Executa o processamento e plotagem
process_and_plot(input_file, output_file, plot_dir, phi_values)

#Essa função vai ler o arquivo para a realizar a simetrização
def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None

    for i in range(17, len(lines)):  # Começar a ler após o cabeçalho de 17 linhas
        line = lines[i].strip()

        if line:
            parts = line.split()

            if len(parts) == 6:
                theta_value = float(parts[3])  # Lendo o valor de theta

            elif len(parts) == 4 and theta_value is not None:
                phi = float(parts[0])
                col1 = float(parts[1])
                col2 = float(parts[2])
                intensity = float(parts[3])
                data.append([phi, col1, col2, theta_value, intensity, True])  # Marcar como original

    df = pd.DataFrame(data, columns=['Phi', 'Col1', 'Col2', 'Theta', 'Intensity', 'IsOriginal'])

    # Debug: Verificar os valores de theta lidos
    print("Valores de Theta lidos:", df['Theta'].unique())

    return df

def save_to_text_file(data_df, intensity_symmetric, output_file_path):
    """
    Salva os dados atualizados em um arquivo de texto no formato solicitado.

    Parameters:
        data_df (DataFrame): O DataFrame com os dados.
        intensity_symmetric (2D array): Intensidades simetrizadas.
        output_file_path (str): O caminho do arquivo de saída.
    """
    with open(output_file_path, 'w') as file:
        # Escrever o cabeçalho
        num_theta = len(data_df['Theta'].unique())
        num_points = len(data_df['Phi'].unique()) * num_theta
        num_phi = len(data_df['Phi'].unique())
        file.write(f"      {num_theta}    {num_points}    0     datakind beginning-row linenumbers\n")
        file.write(f"----------------------------------------------------------------\n")
        file.write(f"MSCD Version 1.00 Yufeng Chen and Michel A Van Hove\n")
        file.write(f"Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720\n")
        file.write(f"Copyright (c) Van Hove Group 1997. All rights reserved\n")
        file.write(f"--------------------------------------------------------------\n")
        file.write(f"angle-resolved photoemission extended fine structure (ARPEFS)\n")
        file.write(f"experimental data for Fe 2p3/2 from Fe on STO(100)  excited with hv=1810eV\n")
        file.write(f"\n")
        file.write(f"provided by Pancotti et al. (LNLS in 9, June 2010)\n")
        file.write(f"   intial angular momentum (l) = 1\n")
        file.write(f"   photon polarization angle (polar,azimuth) = (  30.0,   0.0 ) (deg)\n")
        file.write(f"   sample temperature = 300 K\n")
        file.write(f"\n")
        file.write(f"   photoemission angular scan curves\n")
        file.write(f"     (curve point theta phi weightc weighte//k intensity chiexp)\n")
        file.write(f"      {num_theta}     {num_points}       1       {num_theta}     {num_phi}     {num_points}\n")

        number_of_theta = 0
        for theta in sorted(data_df['Theta'].unique()):
            number_of_theta += 1
            first_row = data_df[data_df['Theta'] == theta].iloc[0]
            file.write(
                f"       {number_of_theta}     {num_phi}       19.5900     {first_row['Theta']:.4f}      1.00000      0.00000\n"
            )
            # Para cada valor de theta, escrever uma linha para cada valor de phi
            theta_data = data_df[data_df['Theta'] == theta]
            sorted_phi_indices = theta_data['Phi'].argsort()  # Ordenar os índices de Phi
            for j, phi in enumerate(theta_data['Phi'].values[sorted_phi_indices]):
                col1 = theta_data.iloc[j]['Col1']
                col2 = theta_data.iloc[j]['Col2']
                intensity = intensity_symmetric[number_of_theta - 1, j]  # Usar o índice correto
                file.write(f"      {phi:.5f}      {col1:.1f}      {col2:.1f}      {intensity:.7f}\n")
                file.write("")  # Linha em branco para separar os blocos de dados

# Leitura dos dados
file_path = '../coeficientes_ajustados.txt'  # Substitua pelo caminho do arquivo

data_df = process_file(file_path)

# Organizar os dados para a simetrização
theta_values = np.sort(data_df['Theta'].unique())  # Garantir que os valores de theta estejam ordenados
phi_values = data_df['Phi'].unique()
n_theta = len(theta_values)
n_phi = len(phi_values)

# Criar uma matriz de intensidades
intensity_values = np.zeros((n_theta, n_phi))
for i, theta in enumerate(theta_values):
    theta_data = data_df[data_df['Theta'] == theta]
    intensity_values[i, :] = theta_data.sort_values(by='Phi')['Intensity'].values

# Aplicar simetrização
intensity_symmetric = fourier_symmetrization(theta_values, phi_values, intensity_values, symmetry)

# Substituir as intensidades originais no DataFrame pelos valores simetrizados
for i, theta in enumerate(theta_values):
    theta_data = data_df[data_df['Theta'] == theta]
    sorted_phi_indices = np.argsort(theta_data['Phi'].values)
    data_df.loc[theta_data.index, 'Intensity'] = intensity_symmetric[i, sorted_phi_indices]

# Salvar os resultados em um arquivo de texto no formato desejado

save_to_text_file(data_df, intensity_symmetric, output_file_path)

def process_file(file_path, sigma=3, rotate_angle=0):
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

    # Suavização da intensidade por Theta
    df['Smoothed_Intensity'] = np.nan
    for theta_value in df['Theta'].unique():
        df_theta = df[df['Theta'] == theta_value].sort_values(by='Phi').copy()
        if len(df_theta) > 1:
            smoothed_intensity = gaussian_filter(df_theta['Intensity'].values, sigma=sigma)
            df.loc[df['Theta'] == theta_value, 'Smoothed_Intensity'] = smoothed_intensity
        elif len(df_theta) == 1:
            df.loc[df['Theta'] == theta_value, 'Smoothed_Intensity'] = df_theta['Intensity'].values[0]

    ####################################################################################################################
    df_plot = df.copy()

    def rotate_phi_for_plot(df_plot, rotation_angle):
        """
        Mesma rotação, mas garante que Phi = 0 esteja presente **apenas** para plotagem.
        Opera diretamente no df_plot de entrada.
        """
        df_plot['Phi'] = (df_plot['Phi'] + rotation_angle) % 360

        if not np.isclose(df_plot['Phi'].min(), 0):
            df_0 = df_plot[df_plot['Phi'] == df_plot['Phi'].min()].copy()
            df_0['Phi'] = 0
            df_plot = pd.concat([df_plot, df_0], ignore_index=True)

        return df_plot
    df_plot = rotate_phi_for_plot(df_plot, rotate_angle)
    # Verificar o intervalo de Phi
    phi_min2 = df_plot['Phi'].min()
    phi_max2 = df_plot['Phi'].max()
    phi_interval2 = phi_max2 - phi_min2

    if phi_interval2 < 360:
        # Encontrar o menor valor de Phi e garantir que ele esteja no intervalo correto
        phi_min2 = df_plot['Phi'].min()

        # Pegar as linhas que têm Phi próximo de phi_min
        df_360 = df_plot[np.isclose(df_plot['Phi'], phi_min2)].copy()
        df_360['Phi'] = 360  # Ajustar para 360°

        df_plot = pd.concat([df_plot, df_360], ignore_index=True)


    if phi_interval2 == 120:
        # Replicação dos dados para cobrir 360 graus
        first_values = df_plot.groupby('Theta').first().reset_index()
        df_plot = df_plot.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df_plot.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original
        # Adicionar os novos valores ao DataFrame
        df_plot = pd.concat([df_plot, last_values], ignore_index=True)

        df_0_120 = df_plot.copy()
        df_0_120['Phi'] = 120 + df_0_120['Phi']
        df_0_120['isOriginal'] = False

        df_240_360 = df_plot.copy()
        df_240_360['Phi'] = 240 + df_240_360['Phi']

        df_plot = pd.concat([df_plot, df_0_120, df_240_360]).reset_index(drop=True)

    if phi_interval2 == 90:
        # Replicação dos dados para cobrir 360 graus
        first_values = df_plot.groupby('Theta').first().reset_index()
        df_plot = df_plot.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df_plot.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original

        # Adicionar os novos valores ao DataFrame
        df_plot = pd.concat([df_plot, last_values], ignore_index=True)
        df_0_90 = df_plot.copy()
        df_0_90['Phi'] = 90 + df_0_90['Phi']

        df_90_180 = df_plot.copy()
        df_90_180['Phi'] = 180 + df_90_180['Phi']

        df_180_270 = pd.concat([df_plot,df_0_90]).reset_index(drop=True)
        df_180_270['Phi'] = 180 + df_180_270['Phi']

        df_plot = pd.concat([df_plot, df_0_90, df_180_270]).reset_index(drop=True)

    if phi_interval2 == 177:
        first_values = df_plot.groupby('Theta').first().reset_index()
        df_plot = df_plot.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original
        df_plot = pd.concat([df_plot, last_values], ignore_index=True)

        df_0_180 = df_plot.copy()
        df_0_180['Phi'] = 180 + df_0_180['Phi']
        df_plot = pd.concat([df_plot, df_0_180]).reset_index(drop=True)

    ####################################################################################################################
    def rotate_phi(df, rotation_angle):
        """
        Rotaciona os valores de Phi no DataFrame garantindo que permaneçam no intervalo [0, 360).

        Parâmetros:
        - df (pd.DataFrame): DataFrame com a coluna 'Phi'.
        - rotation_angle (float): Ângulo de rotação em graus.

        Retorna:
        - pd.DataFrame: DataFrame com Phi rotacionado.
        """
        df['Phi'] = (df['Phi'] + rotation_angle) % 360
        return df
        # Verificar o intervalo de Phi

    df = rotate_phi(df, rotate_angle)

    phi_min = df['Phi'].min()
    phi_max = df['Phi'].max()
    phi_interval = phi_max - phi_min

    if phi_interval < 360:
        # Encontrar o menor valor de Phi e garantir que ele esteja no intervalo correto
        phi_min = df['Phi'].min()

        # Pegar as linhas que têm Phi próximo de phi_min
        df_360 = df[np.isclose(df['Phi'], phi_min)].copy()
        df_360['Phi'] = 360  # Ajustar para 360°

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

    if phi_interval == 177:
        first_values = df.groupby('Theta').first().reset_index()
        df = df.groupby('Theta', group_keys=False).apply(lambda x: x.drop(x.index[0]))
        last_values = df.groupby('Theta').last().reset_index()  # Pegar os últimos valores
        last_values['Phi'] = first_values['Phi']  # Substituir pelo valor do primeiro Phi original
        df = pd.concat([df, last_values], ignore_index=True)

        df_0_180 = df.copy()
        df_0_180['Phi'] = 180 + df_0_180['Phi']
        df = pd.concat([df, df_0_180]).reset_index(drop=True)

    return df, df_plot


# Função para interpolar os dados
def interpolate_data(df, resolution=1000):
    # Definir uma grade regular para a interpolação
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['Smoothed_Intensity'] # Usar a intensidade suavizada

    # Criando uma grade de pontos onde queremos interpolar
    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    # Criando uma grade de malha para interpolação
    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    # Realizar a interpolação
    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')

    return phi_grid, theta_grid, intensity_grid

from scipy.ndimage import gaussian_filter
# Função para gerar o gráfico polar
def plot_polar_interpolated(df, resolution=500):

    # Interpolar os dados
    plt.ion()
    matplotlib.use('TkAgg')
    phi_grid, theta_grid, intensity_grid = interpolate_data(df, resolution)

    #intensity_grid = gaussian_filter(intensity_grid, sigma=3) # Você pode suavizar a grade interpolada também
    # Criando o gráfico polar
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10,8), dpi=100)



    norm = mcolors.PowerNorm(gamma=1.0, vmin=intensity_grid.min(), vmax=intensity_grid.max())

    # Plotando a intensidade interpolada
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot', norm=norm)

    # Definir o limite máximo do eixo theta com base no maior valor de theta nos dados
    max_theta = df['Theta'].max()
    ax.set_ylim(0, np.radians(max_theta))

    # Adiciona rótulos para os ângulos theta
    theta_ticks = np.linspace(0, max_theta, num=6)
    ax.set_yticks(np.radians(theta_ticks))
    ax.set_yticklabels([f'{int(tick)}°' for tick in theta_ticks])

    ax.set_xlabel('', fontsize=36)
    ax.set_ylabel('', fontsize=36, labelpad=1, rotation=360)
    ax.yaxis.set_label_coords(0.8, 0.93)

    # Adicionando a barra de cores
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label('', fontsize=36)

    plt.xticks(fontsize=11)
    plt.yticks(fontsize=0)
    plt.draw()
    plt.pause(600)
    plt.close()


# Caminho do arquivo de dados
file_path = '../simetrizados.txt'





def save_to_txt_with_blocks(df, file_name):
    """
    Salva os dados organizados em blocos, onde cada bloco corresponde a um valor de θ (Theta),
    removendo o valor de Phi = 360 caso sua intensidade seja igual à do primeiro valor de Phi do bloco.
    """
    # Filtrar apenas os dados originais
    df_original = df[df['IsOriginal']]

    num_theta = df_original['Theta'].nunique()
    num_phi = df_original['Phi'].nunique() - 1
    num_points = num_theta*num_phi

    with open(file_name, 'w',newline="\n") as file:
        # Cabeçalho inicial
        file.write(f"   323    17    0     datakind beginning-row linenumbers\n")
        file.write(f"----------------------------------------------------------------\n")
        file.write(f"MSCD Version 1.00 Yufeng Chen and Michel A Van Hove\n")
        file.write(f"Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720\n")
        file.write(f"Copyright (c) Van Hove Group 1997. All rights reserved\n")
        file.write(f"--------------------------------------------------------------\n")
        file.write(f" angle-resolved photoemission extended fine structure (ARPEFS)\n")
        file.write(f" experimental data for Pd 3d3/2 from W(100)  excited with hv=1810eV\n")
        file.write("\n")
        file.write(f" provided by Pancotti et al. (LNLS in 12, June 2007)\n")
        file.write(f"   intial angular momentum (l) = 2\n")
        file.write(f"   photon polarization angle (polar,azimuth) = (  30.0,   0.0 ) (deg)\n")
        file.write(f"   sample temperature = 300 K\n")
        file.write("\n")
        file.write(f"   photoemission angular scan curves\n")
        file.write(f"     (curve point theta phi weightc weighte//k intensity chiexp)\n")
        file.write(f"      {num_theta}    {num_points}       1      {num_theta}     {num_phi}    {num_points}\n")

        # Número do bloco de θ
        number_of_theta = 0

        for theta in sorted(df_original['Theta'].unique()):
            number_of_theta += 1
            subset = df_original[df_original['Theta'] == theta].sort_values(by='Phi')
            first_intensity = subset.iloc[0]['Intensity']
            last_intensity = subset.iloc[-1]['Intensity']
            if number_of_theta < 10:
                theta_format = f"       {number_of_theta}"
            else:
                theta_format = f"      {number_of_theta}"
            file.write(
                f"{theta_format}     {num_phi}      9.85000      {theta:1.0f}      1.00000      0.00000\n")

            for i, row in enumerate(subset.itertuples(index=False)):
                if not (i == len(subset) - 1 and row.Intensity == first_intensity):
                    if row.Phi < 10:
                        phi_format = f"{row.Phi:7.5f}"  # Exemplo: " 0.12345"
                    elif row.Phi < 100:
                        phi_format = f"{row.Phi:7.4f}"  # Exemplo: "12.3456"
                    else:
                        phi_format = f"{row.Phi:7.3f}"  # Exemplo: "123.456"
                    if row.Col1 < 100000:  # Menos de 100000
                        col1_format = f"{row.Col1:1.1f}"  # Exemplo: " 98755.2"
                    else:  # Acima de 100000
                        col1_format = f"{row.Col1:1.0f}."  # Exemplo: "100697.0" (com ponto
                    if row.Col2 < 100000:  # Menos de 100000
                        col2_format = f"{row.Col2:1.1f}"  # Exemplo: " 98755.2"
                    else:  # Acima de 100000
                        col2_format = f"{row.Col2:1.0f}."  # Exemplo: "100697.0" (com ponto
                    file.write(f"      {phi_format}      {col1_format}      {col2_format}   {row.Intensity: 10.7f}\n")


# Processar os dados
df, df_plot = process_file(file_path, 1, rotate_angle)


plot_polar_interpolated(df_plot)

save_to_txt_with_blocks(df, arquivo_saida)


