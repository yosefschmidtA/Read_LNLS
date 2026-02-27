#!/bin/bash

cp OB.out  tempA.txt;
cp OA.out  tempB.txt;

# Construindo os Vetores Necessarios

coluna1=($(awk '{print $1}' tempA.txt | xargs)); 
coluna3=($(awk '{print $3}' tempA.txt | xargs)); 
coluna4=($(awk '{print $4}' tempA.txt | xargs)); 
coluna5=($(awk '{print $5}' tempA.txt | xargs));

coluna2A=($(awk '{print $2}' tempA.txt | xargs));
coluna2B=($(awk '{print $2}' tempB.txt | xargs));

linhas=$(cat tempA.txt | wc -l);  # conta o numero de linhas  no arquivo 
linhas=$(echo "$linhas -1" | bc); # subtrai um, pois o vetor comeca no indice 0

# Faz a mistura ponderada para os varios modelos

for i in `seq 0 $linhas`;
       do
             coluna2[$i]=$(echo "print ${coluna2A[$i]}*0.78+${coluna2B[$i]}*0.22" | python);
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

        cat > mistura.out << acabou
   223    26     0     datakind beginning-row multi-curves
-----------------------------------------------------------------
Calculation of photoelectron diffraction and dichroism
MSCD Parallel  Version 1.37 Yufeng Chen and Michel A Van Hove
Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720
Copyright (c) Van Hove Group 1997-1998. All rights reserved
-----------------------------------------------------------------
 angle-resolved photoemission extended fine structure (ARPEFS)
 multiple scattering calculation of Ga2O3(100)
 calculated by Yosef (UFJ-PPGCAS) on Apr 27, 2025

   initial angular momentum (li) = 0   msorder=  6   raorder= 4
   photon polarization angle (polar,azimuth) =(   0.0   0.0 ) (deg)

   radius, depth and lattice constant=   9.0,  11.3 and   3.04 angstrom
   cluster size= 122 atoms and spacings=  1.43  1.82  1.70  1.55 angs
   inner potential=  0.0 V  debye and sample temperature= 105, 300 K
   number of valence electrons=  24     bandgap energy=   4.70 eV
   density of bulk=   5.88 g/cm3     molecular weight= 187.44 amu
   effective weight for kind 1-3 =   16.0   69.7   55.8 amu
   half aperture angle=    3.0 deg            pathcut=     0.0010

   photoemission azimuthal scan curves
     parameters: curve point theta phi weightc weighte
     columns: phi intensity background chical chiexp
    20  2400    1   20  120 2400  ncurve npoint nk ntheta nphi nangle
acabou

cat temp.txt >> mistura.out ;

rm -rf temp.txt;
rm -rf tempA.txt;
rm -rf tempB.txt;

echo "Ok!"
