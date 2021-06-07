import MDAnalysis
import os
from os import listdir
from os.path import isfile, join
import mdtraj as md
import math
import numpy as np
mypath='.'


def loadstructures():
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-4:]=='.pdb' and '_ali' in f]
    atoms={}
    for f in onlyfiles :
        print(f)
        atoms.update({f[:-4]:md.load(f)})
    return atoms

def Calpha_coordinate(traj):
    try :return traj.atom_slice(traj.topology.select('name CB')  ).xyz[0][0]
    except :
        return False

def is_aligned (a, b , c):
    if 160  <= coordangle (a,b,c)<= 200 :
        return True
    else : return False

def not_overlap (a, b ):

    if  coordistance (a,b)> 0.45 :
        return True
    else : return False

def not_clashing (traja, trajb ):
    for a in traja.atom_slice(traja.topology.select('all')  ).xyz[0] :
        for b in  trajb.atom_slice(trajb.topology.select('all')  ).xyz[0] :
            if  coordistance (a,b)< 0.1 :
                return False
    return True
def close_enough (a, b ):

    if  0.45 < coordistance(a,b) <0.65 :
        return True
    else : return False


def coordistance (a,b) :
    return math.sqrt((a[0]- b[0])**2 + (a[1]- b[1])**2 + (a[2]- b[2])**2  )

def coordangle (a,b,c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)



to_keep3=[]
to_keep4=[]
sidechaindict=loadstructures()
for i , (sidechain1, traj1) in  enumerate(sidechaindict.items()) :
    print(i , len(to_keep3),  len(to_keep4))
    coord1= Calpha_coordinate(traj1)
    for j , (sidechain2, traj2) in  enumerate(sidechaindict.items())  :
        coord2= Calpha_coordinate(traj2)
        if isinstance(coord1,(list, tuple, np.ndarray)) and isinstance(coord2,(list, tuple, np.ndarray)) and sidechain1!=sidechain2 and not_clashing (traj1, traj2 ) and close_enough (coord1, coord2 ) :
            for k , (sidechain3,traj3 )  in  enumerate(sidechaindict.items()) :
                coord3=Calpha_coordinate(traj3)
                if  isinstance(coord3,(list, tuple, np.ndarray)) and sidechain2!=sidechain3 and  not_clashing (traj2, traj3 )  and close_enough (coord2, coord3 ) and is_aligned (coord1, coord2 , coord3) :
                    to_keep3.append([sidechain1,sidechain2,sidechain3])
                    print( [sidechain1,sidechain2,sidechain3])
                    #if sidechain1[:4]==sidechain2[:4] and sidechain1[:4]==sidechain3[:4] :
                        #os.system('echo HEADER %s %s %s  >sameantibodychain%s.pdb ' %(sidechain1,sidechain2,sidechain3,i) )
                        #for file in [sidechain1,sidechain2,sidechain3] :
                        	#os.system('cat %s.pdb >> sameantibodychain%s.pdb' %(file,i) )
                    #else :
                        #os.system('echo HEADER %s %s %s  > out%s.pdb' %(sidechain1,sidechain2,sidechain3,i) )
                        #for file in [sidechain1,sidechain2,sidechain3] :
                        	#os.system('cat %s.pdb >> out%s.pdb' %(file,i) )
                    for l , (sidechain4,traj4 )  in   enumerate(sidechaindict.items()) :
                        coord4=Calpha_coordinate(traj4)
                        if isinstance(coord4,(list, tuple, np.ndarray)) and sidechain3!=sidechain4  and  not_clashing (traj3, traj4 )  and close_enough (coord3, coord4 ) and is_aligned ( coord2 , coord3, coord4):
                            to_keep4.append([sidechain1,sidechain2,sidechain3,sidechain4])
                            #os.system('echo HEADER %s %s %s %s  > 4out%s.pdb' %(sidechain1,sidechain2,sidechain3,sidechain4,i) )
                            #for file in [sidechain1,sidechain2,sidechain3,sidechain4] :
                                #os.system('cat %s.pdb >> 4out%s.pdb' %(file,i) )

print(to_keep3)
print(to_keep4)


for i, combinaison in enumerate(to_keep3) :
        os.system('rm out%s.pdb' %(i) )
        sidechain1= combinaison[0]
        sidechain2= combinaison[1]
        sidechain3= combinaison[2]
        if sidechain1[:4]==sidechain2[:4] and sidechain1[:4]==sidechain3[:4] :
            outfile = 'sameantibodychain'
        else :     outfile = 'out'
        os.system('echo HEADER %s %s %s  > %s%s.pdb' %(combinaison[0] ,combinaison[1] ,combinaison[2] ,outfile,i) )
        for file in  combinaison :
            os.system('cat %s.pdb >> %s%s.pdb' %(file,outfile, i) )



for j, combinaison in enumerate(to_keep4) :
        sidechain1= combinaison[0]
        sidechain2= combinaison[1]
        sidechain3= combinaison[2]
        sidechain4= combinaison[3]
        os.system('rm out%s.pdb' %(i) )

        if sidechain1[:4]==sidechain2[:4] and sidechain1[:4]==sidechain3[:4] and sidechain1[:4]==sidechain4[:4] :
            outfile = '4sameantibodychain'
        else :     outfile = '4out'
        os.system('echo HEADER %s %s %s %s  > %s%s.pdb' %(combinaison[0] ,combinaison[1] ,combinaison[2] ,combinaison[3] ,outfile,i) )
        for file in  combinaison :

            os.system('cat %s.pdb >> %s%s.pdb' %(file,outfile, i) )
