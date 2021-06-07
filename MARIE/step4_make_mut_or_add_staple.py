
from library_scanning_ubc import mutate_residues
from library_scanning_ubc import add_staple
from library_scanning_ubc import add_fluo
from library_scanning_ubc import prep_filexfold
import MDAnalysis as mda
from os import listdir
from os.path import isfile, join
import math
import sys

#input=sys.argv[1]

def loadinput(mypath='.'):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-4:]=='.pdb' and 'output' in f and 'pep' not in f]
    return onlyfiles



def coordistance (a,b) :
    return math.sqrt((a[0]- b[0])**2 + (a[1]- b[1])**2 + (a[2]- b[2])**2  )
def loadstructures(mypath='../'):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-4:]=='.pdb' and '_ali' in f]
    atoms={}
    for f in onlyfiles :
        #print(f)

        atoms.update({f[:-4]:mda.Universe(mypath+f)})
    return atoms

def solventexposedala (structure1,structure2 ) :
    u2 = mda.Universe(structure2).atoms.positions
    notinteracting=[]
    for res in  mda.Universe(structure1).select_atoms('resname ALA').chain('A').residues :
            exposed=True
            for atom_pos2 in u2 :
                if coordistance (res.center_of_mass(), atom_pos2) < 15 :
                    exposed=False
            if exposed==True :
                notinteracting.append(res.resid)
    return notinteracting




for fi in  loadinput() :
    print('dealing with file :' + fi)
    prep_filexfold(fi, ['A'] )
    filename=fi[:-4] + '_prep.pdb'

    final= mda.Universe(filename)
    mutationspoints=final.select_atoms('resname ALA and name CB').positions
    resindex = final.select_atoms('resname ALA and name CB').resids
    print(resindex)
    sidechains=  loadstructures()
    print()
    pos_mutations =[]
    for j, pos in enumerate(mutationspoints) :
        print(pos)
        for i , (sidechain, traj) in  enumerate(sidechains.items()) :
            new_pos= traj.select_atoms('name CB').positions
            if len(new_pos)>0 :
                new_pos=new_pos[0]
                #print(coordistance ( pos, new_pos))
                if coordistance ( pos, new_pos) <3 :
                    print('mutation possible'  + str(resindex[j]) +  traj.select_atoms('name CB').resnames[0])
                    if ['P' , resindex[j] ,  traj.select_atoms('name CB').resnames[0]] not in pos_mutations :
                        pos_mutations.append(['P' , resindex[j] ,  traj.select_atoms('name CB').resnames[0]] )


    for i , mut in enumerate(pos_mutations) :
        print(mut,filename )

        #mutate_residues( 'output.pdb', [mut],  '%s-%s' % (input, i) )
        #add_staple('%s-%s_full.pdb' % (input, i) ,chain='P')
        print( '%s-%s' % (filename[:-4], i))

        mutate_residues( filename , [mut],  '%s-%s' % (filename[:-4], i) )
        # .add_staple('%s-%s_full.pdb' % (input, i) ,chain='P')



def oneresscanallbutone(filename ,res='ALA', fileout=False):
    if fileout==False :
        fileout=filename[:-4]
    for j in sequence :
        mutations = [ ['D' , k , 'ALA']  for k in sequence if k!=j ]
        mutate_residues( filename ,mutations, '%s-scan-%s-%s' % (res , j , fileout ))


def oneresscan(filename ,res='ALA', fileout=False ):
    if fileout==False :
        fileout=filename[:-4]
    for j in surface_facing_res :
        mutations = [['D' , j , res]]
        mutate_residues( filename,mutations,  '%s-scan-%s-%s' % (res , j , fileout) )
        add_staple('%s-scan-%s-%s_full.pdb' % (res , j , fileout), chain='P')


filename='output.pdb'
