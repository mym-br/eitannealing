#!/bin/bash
for i in {1..16}
do
   echo $i
   for j in 0.875 0.5 0.05 0.025
      do
       echo $j;
       build/Release/eitannealingtest circular_2D.msh ueta_circuito_basal_0XXX_CORRENTES_10mA.txt ueta_circuito_0XXX_AMPLITUDES_10mA.txt $j 
      done
done
