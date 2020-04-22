pdb2gmx -f shorthelix.pdb -o shelix.gro -p shelix.top -i shelix.posre.itp -ignh -ss
editconf -f shelix.gro -o shelix.box.gro -d 1.2 -bt octahedron
genbox -cp shelix.box.gro -cs spc216.gro -o shelix.solv.gro -p shelix.top
grompp -f ions.mdp -c shelix.solv.gro  -p shelix.top -o ions.tpr
genion -s ions.tpr -o shelix.ready.gro -p shelix.top -neutral
grompp -f ions.mdp -c shelix.ready.gro -p shelix.top -o shelix.min.tpr
mdrun -deffnm shelix.min -v
grompp -f nvt.mdp -c shelix.min.gro -p shelix.top -o shelix.nvt.tpr
mdrun -deffnm shelix.nvt 
grompp -f npt.mdp -c shelix.nvt.gro -p shelix.top -o shelix.npt.tpr
mdrun -deffnm shelix.npt 
