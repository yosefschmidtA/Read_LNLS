import numpy as np
import pandas as pd
import argparse


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
    lista_temporaria = []

    for i in range(26, len(lines)):
        line = lines[i].strip()

        if "fitted parameters (" in line:
            break

        if line:
            parts = line.split()

            if len(parts) == 7:
                # Encontrou um novo Theta, processa os dados acumulados
                if lista_temporaria:
                    intensi2_values = [item[1] for item in lista_temporaria]
                    media_intensi2 = sum(intensi2_values) / len(intensi2_values)

                    for phi, intensi1, intensi2, intensityexp in lista_temporaria:
                        intensitycal = (intensi1 - media_intensi2) / media_intensi2
                        data.append([phi, intensitycal, theta_value, intensityexp, True])

                    lista_temporaria = []

                theta_value = float(parts[3])

            elif len(parts) == 5 and theta_value is not None:
                phi = float(parts[0])
                intensi1 = float(parts[1])
                intensi2 = float(parts[2])
                intensityexp = float(parts[4])
                lista_temporaria.append((phi, intensi1, intensi2, intensityexp))

    # Depois do loop, processa o último conjunto também
    if lista_temporaria:
        intensi2_values = [item[1] for item in lista_temporaria]
        media_intensi2 = sum(intensi2_values) / len(intensi2_values)

        for phi, intensi1, intensi2, intensityexp in lista_temporaria:
            intensitycal = (intensi1 - media_intensi2) / media_intensi2
            data.append([phi, intensitycal, theta_value, intensityexp, True])

    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])

    # Calcula R-factor
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)

    for theta, r_factor in zip(angulos, resultados):
        print(f'Theta: {theta}, R-factor: {r_factor}')

    print(f'R-factor médio: {r_factor_medio}')
    print(f'R-factor total: {r_factor_total}')
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

    return r_factor_total  # Retorna r_factor_total

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calcular R-factor a partir de um arquivo de saída.")
    parser.add_argument("arquivo", type=str, help="Arquivo de entrada")
    parser.add_argument("valor1", type=str, help="Valor 1")
    parser.add_argument("valor2", type=str, help="Valor 2")
    parser.add_argument("saida", type=str, help="Arquivo de saída")  # Adicionando argumento de saída

    args = parser.parse_args()

    r_factor_total = process_file(args.arquivo)

    # Escrever no arquivo de saída
    with open(args.saida, "a") as f:
        f.write(f"{args.valor1}  {args.valor2}  {r_factor_total}\n")



