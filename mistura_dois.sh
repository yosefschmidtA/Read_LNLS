#!/bin/bash

export LC_NUMERIC="C"
# Verificando se o passo foi passado como argumento
if [ $# -ne 1 ]; then
    echo "Uso: $0 <passo>"
    exit 1
fi

passo=$1

cp OB.out tempA.txt;
cp OA.out  tempB.txt;

# Construindo os Vetores Necessários

coluna1=(`awk '{print $1}' tempA.txt | xargs`); 
coluna3=(`awk '{print $3}' tempA.txt | xargs`); 
coluna4=(`awk '{print $4}' tempA.txt | xargs`); 
coluna5=(`awk '{print $5}' tempA.txt | xargs`);

coluna2A=(`awk '{print $2}' tempA.txt | xargs`); 
coluna2B=(`awk '{print $2}' tempB.txt | xargs`);

linhas=`cat tempA.txt | wc -l`;  # conta o numero de linhas no arquivo 
linhas=`echo "$linhas -1" | bc`; # subtrai um, pois o vetor comeca no indice 0

# Limpando arquivo de saída anterior
rm -f saida_mistura.s1

# Loop para testar combinações de fator1 e fator2
for fator1 in $(seq 0.00 $passo 1.00); do
    fator2=$(echo "1 - $fator1" | bc)  # fator2 é sempre 1 - fator1

    # Faz a mistura ponderada para os varios modelos
    for i in `seq 0 $linhas`
    do
        coluna2[$i]=`echo "print ${coluna2A[$i]}*${fator1}+${coluna2B[$i]}*${fator2}" | python`;
        echo "Construindo linha $i de $linhas";
    done;

    for i in `seq 0 $linhas`;
    do
        case "${coluna2A[$i]}" in
        120) 
            echo "  ${coluna1[$i]}   ${coluna2[$i]}   ${coluna3[$i]}   ${coluna4[$i]}   ${coluna5[$i]}   0.0  ------------------" >> temp.txt;;      
        *)
            echo "  ${coluna1[$i]}   ${coluna2[$i]}   ${coluna3[$i]}   ${coluna4[$i]}   ${coluna5[$i]}" >> temp.txt;;
        esac 
    done;

    # Gerando o arquivo de saída para o cálculo do R
    cat > saida_mistura.s1 << acabou
   2400    20     0     datakind beginning-row multi-curves
-----------------------------------------------------------------
Calculation of photoelectron diffraction and dichroism
MSCD Parallel  Version 1.37 Yufeng Chen and Michel A Van Hove
Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720
Copyright (c) Van Hove Group 1997-1998. All rights reserved
-----------------------------------------------------------------
 angle-resolved photoemission extended fine structure (ARPEFS)
 multiple scattering calculation of SiC(0001)-Grafeno
 calculated by LHLima on Jul 21, 2011

   initial angular momentum (li) = 1   msorder=  8   raorder= 4
   photon polarization angle (polar,azimuth) =(  30.0   0.0 ) (deg)

   radius, depth and lattice constant=   9.0,  16.0 and   3.08 angstrom
   cluster size= 329 atoms and spacings=  0.00  0.10  1.96  0.55 angs
   inner potential= 11.7 V  debye and sample temperature= 625, 300 K
   electron wave attenuation due to inelastic process not considered
   density of bulk=   3.21 g/cm3     molecular weight=  40.10 amu
   effective weight for kind 1-3 =   12.0   28.1    0.0 amu
   half aperture angle=    4.0 deg            pathcut=     0.0010

   photoemission azimuthal scan curves
     parameters: curve point theta phi weightc weighte
     columns: phi intensity background chical chiexp
    20  2400    1   20  120 2400  ncurve npoint nk ntheta nphi nangle
acabou

    cat temp.txt >> saida_mistura.s1;

    # Executando o cálculo do R após a mistura
    python3 calcularRMistura.py saida_mistura.s1 $fator1 $fator2 resultado_mistura.txt

    # Limpando arquivos temporários
    rm -rf temp.txt;
    rm -rf tempA.txt;
    rm -rf tempB.txt;

    echo "Cálculo do R finalizado para fator1=$fator1 e fator2=$fator2"
done

echo "Processo finalizado!"

