#!/bin/bash
####### BEGIN-CHANGEME
source /usr/local/gromacs/bin/GMXRC
#export GMXLIB=/home/marie/cov2-s-pep/MARIE/FFfiles/
####### END-CHANGEME
gpu=1
echo "Using GPU id 1. Change the variable gpu to a diffrerent id if required"

####Minimisation
rm \#*
gmx -quiet grompp -f ions.mdp -c  shelix.ready.gro -p  shelix.top -maxwarn 1  -o  em.tpr >& grompp.trash 
gmx -quiet  mdrun -gpu_id $gpu -v -c shelix.minimised.gro -s em.tpr >& minimisation.output &
echo "GROMACS is running the minimisation step. This should not take long..."
wait
echo "Done! Standard output has been redirected to minimisation.output"
rm \#*
sleep 1
###equilibration 1 NVT
gmx -quiet grompp -f nvt.mdp -c shelix.minimised.gro -p shelix.top -o eq.tpr  -r shelix.minimised.gro -maxwarn 1 >& grompp.trash 
gmx -quiet mdrun -gpu_id $gpu -s eq.tpr -v -nt 2 -nice 0 -c shelix.nvt.gro  -ntmpi 1 >& nvt.output &
echo "GROMACS is running a short NVT equilibration.."
sleep 2
echo "Boy, this is taking longer... While we wait, let's meet some dutch painters!"
sleep 7
echo "Sorry, I did not find any dutch painter of relevance other than a Van Gue? Van Gal? Van der Waals? Let me check something..."
wait
echo "Done! Standard output has been redirected to nvt.output"
rm \#*
sleep 1
gmx -quiet grompp -f npt.mdp -c shelix.nvt.gro -p shelix.top -o eq.tpr -r shelix.nvt.gro -maxwarn 1 >& grompp.trash
gmx -quiet mdrun -gpu_id $gpu -s eq.tpr -v -nt 2 -nice 0 -c shelix.npt.gro  -ntmpi 1 >& npt.output &
echo "GROMACS is running a short NPT equilibration.."
sleep 2
echo "Where was I? Ah yes, no dutch painters! Let's check the spanish wikipedia page for flemish painting instead!"
sleep 2                                                  
echo "---Late Gothic and Early Renaissance---" 
sleep 1
echo "Jan van Eyck (1390-1441)"
sleep 1
echo "Hans Memling"
sleep 1
echo "Rogier van der Weyden (1399-1464)"
sleep 1
echo "Hugo van der Goes (1440-1482)"
sleep 1
echo "Robert Campin (1378-1444)"
sleep 1
echo "Jheronimus Bosch, llamado en España El Bosco (1450-1516)"
sleep 1
echo "Quentin Massys (m. 1530)"
sleep 1
echo "Joachim Patinir (1480-1524)"
sleep 1
echo "Antonio Moro (m. 1576)"
echo "Wait, there are even more!"
sleep 3
echo "---Baroque---"
echo "Rubens (1577-1640)"
sleep 1
echo "Anton van Dyck (1599-1641)"
sleep 1
echo "Jacob Jordaens (1593-1678)"
sleep 1
echo "Frans Hals (1584-1666)"
sleep 1
echo "Rembrandt (1607-1669)"
sleep 1
echo "Vermeer de Delft (1632-1675)"
sleep 1
echo "Jacob Ruysdael (1628-1682)"
sleep 3
echo "That is something else, isn't it?"
wait
echo "Ah yes! I almost forgot that the GROMACS NPT thing is done! Standard output has been redirected to npt.output"
rm \#*                               
sleep 1
gmx -quiet grompp -f md.mdp -c shelix.npt.gro -p shelix.top -o md.tpr -r shelix.npt.gro -maxwarn 1 >& grompp.trash
gmx -quiet mdrun -gpu_id $gpu -s md.tpr -v -nt 2 -nice 0 -c shelix.md.gro  -ntmpi 1 >& md.output &
echo "Now we are talking business! GROMACS is runing 1 ns MD simulation"
wait
natoms=`head shelix.minimised.gro -n 2 | tail -n 1`
nsday=`grep 'Performance' md.output | awk '{print $2}'`

echo "GROMACS has finished runing. Standard output has been redirected to md.output"
echo "Final Performance:"
echo "n atoms: $natoms prod: $nsday ns/day"
echo "-----WARNING-----"
echo "Something is wrong with the peripheral connecting the keyboard to the chair!"
echo "Check the results carefully!" 
echo "Happy Benchmarking!"

