#/bin/bash

parm7=input_nosolty.prm7
name=`basename $parm7 .prm7`
echo $parm7 
echo $name

$AMBERHOME/bin/ante-MMPBSA.py -p $parm7 -c $name.comp.prm7 -r $name.rec.prm7 -l $name.lig.prm7 -s :WAT,Cl-,Na+ -n :203-216
$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa.in -o ${name}_mmpbsa.dat -sp $parm7 -cp $name.comp.prm7 -rp $name.rec.prm7 -lp $name.lig.prm7 -y $name.nc
