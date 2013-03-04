#!/bin/bash
for i in {1..16}
do
   echo $i
   for j in 0.0001 0.00002 0.000005 0.0000001  
      do
       echo $j;
       build/Release/eitannealingtest circular_2D.msh ueta_circuito_basal_0XXX_CORRENTES_10mA.txt ueta_circuito_0XXX_AMPLITUDES_10mA.txt $j y
      done
done
