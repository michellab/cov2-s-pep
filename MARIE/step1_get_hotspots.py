from os import listdir
from os.path import isfile, join
import MDAnalysis
import mdtraj
import numpy as np
import os
import matplotlib.pylab as plt
mypath='backup/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-4:]=='.pdb']



print(onlyfiles)



chain_RBD=resid='resid 330-528'
chain_antibody='resid 1-250'
def get_interacting_residues_selections(input ) :
    u = MDAnalysis.Universe(input)

    interactingantibody = u.select_atoms('%s and around 4 %s'%(  chain_antibody , chain_RBD))
    antiresid=[]
    for at in interactingantibody.atoms :
        if at.resid not in antiresid :
            antiresid.append(at.resid)

    #interactingantibody.write( '%s-intleractions.pdb' %( input[:-4] ))

    interactingRBD = u.select_atoms(' %s and around 4 %s'%( chain_RBD, chain_antibody ))
    RBDresid=[]
    for at in interactingRBD.atoms :
        if at.resid not in RBDresid :
            RBDresid.append(at.resid)
    '''
    for index in RBDresid :
        u.select_atoms('%s and around 4 resid %s'%(  chain_antibody ,index ))
        if at.resid not in antibodysidechain:
            RBDresid.append(at.resid, at.resname)
    '''
    #print(RBDresid)
    return RBDresid, antiresid

def make_hotmap(plain_RBD,complex) :
    RBDresidall=[]
    uRBDalone = MDAnalysis.Universe(plain_RBD)
    rbd=mdtraj.load_pdb(plain_RBD)
    n_atoms=rbd.n_atoms
    bfactors=np.ones(n_atoms)
    for input in complex:
        print(input)
        RBDresid, antiresid= get_interacting_residues_selections(mypath + input )
        RBDresidall.extend(RBDresid)

        for resid in RBDresid :
            print(input ,resid)
            print(rbd.topology.select('resSeq %s' %resid ))
            for atom in rbd.topology.select('resSeq %s' %resid ):
                bfactors[atom]+=1
                print(input ,resid ,  atom  )


    d = {item:RBDresidall.count(item) for item in RBDresidall }
    print(d)
    lists = sorted(d.items()) # sorted by key, return a list of tuples

    x, y = zip(*lists) # unpack a list of pairs into two tuples

    plt.plot(x, y)
    plt.savefig('out.png')

    #print( bfactors)
    rbd.save_pdb('out.pdb', force_overwrite=True, bfactors=bfactors)
    antiresid=[]
    out=[]
    for key in d.keys() :
        if  d.get(key) > 3:

            for input in complex:



                u = MDAnalysis.Universe(mypath +input)
                chainids = list(dict.fromkeys( at.chainID for at in u.select_atoms('all')))
                print(chainids)
                interactingantibody = u.select_atoms('%s and around 4 resid %s and protein'%(  chain_antibody , key))
                antiresid.clear()

                for at in interactingantibody.atoms :
                    print(at.resid)
                    if at.resid not in antiresid :
                        antiresid.append(at.resid)
                        out.append([key , input, at.resname])
                        print([key , input, at.resname])
                        #to_pdb=u.select_atoms('resid %s  and around 4 resid %s'%(  at.resid ,  key ))
                        #to_pdb=u.select_atoms(' same resid  as (resid %s  and around 4 resid %s) '%(  at.resid ))
                        print( 'resid %s  '%(  at.resid  ))
                        to_pdb=u.select_atoms('resid %s  '%(  at.resid  ))

                        #to_pdb.write( '%s-%s-%s.pdb' %( key, input[:-4],at.resid  ))
                        print(at.chainID )
                        chain = chr(chainids.index(at.chainID ) +ord('A'))
                        print(chain )
                        #to_pdb=u.select_atoms(' same resid  as (resid %s  and around 4 resid %s) '%(   at.resid ,  key ) )
                        #to_pdb.write( '%s-%s-%s.pdb' %( key, input[:-4],at.resid  ))
                        to_pdb.write( 'temp-%s-%s-%s.pdb' %(key, input[:-4],at.resid  ) )
                        os.system('grep \ %s\ \   temp-%s-%s-%s.pdb > %s-%s-%s.pdb'  %(chain  ,key, input[:-4],at.resid , key, input[:-4],at.resid  )   )

    #print(out)

make_hotmap('box.pdb',onlyfiles )
