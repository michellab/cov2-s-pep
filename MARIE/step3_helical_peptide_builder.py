import numpy as np
import subprocess
import copy
import mdtraj
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mda
import os
from os import listdir
from os.path import isfile, join

import sys
import math
#rmsd to the initial -> The bigger the more structures are accepted
#rmsd_cutoff=[1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0]
rmsd_cutoff=[4.5,5.0,5.5,6.0]
# Evaluation of the clashes --> The bigger the 'worst' and so more structures are accepted
clashes_cutoff=80
#sequences start from sequence 'start' to sequence 'end'
start=0
end=174
#end=len(sequences)

class Flat_File :


    def __init__(self, id = str(), path = '.') :
        self.id = id
        self.path = path
        self.lines = list()


    def read_file(self, path = '.') :
        ''' Read a flat file. Assigns a lines list to a lines attribute. This
        fonction is used by CSV and PDB files.
        '''
        if path != '.' :
            self.path = path
        f = open(self.path, "r")
        lines = f.readlines()
        f.close()
        self.lines = lines

    def split_PDBfile_by_chains(self, output_dir = '.', chains = 'all', all_sections = True ) :
        ''' Split a pdb file in different pdb files by chains. data is a list of
        pdb file lines. chains must be a list of PDB ids (e.g. ['A', 'B'])
        '''
        pdblines = self.lines
        # file split :
        initial_sections = list()
        dict_chains = dict()
        final_sections = list()
        i = 0
        while i < len(pdblines) :
            line = pdblines[i]
            if line[0:4] != 'ATOM' and line[0:3] != 'TER' :
                initial_sections.append(line)
                i += 1
            else :
                break

        while i < len(pdblines) :
            line = pdblines[i]
            possible_sections = ['ATOM  ', 'ANISOU', 'TER   ', 'HETATM']
            if line[0:6]in possible_sections:
                chain_id = line[21]
                if not(chain_id in dict_chains) :
                    dict_chains[chain_id] = [line]
                else :
                    dict_chains[chain_id].append(line)
                i += 1
            else :
                i += 1
        while i < len(pdblines) :
            line = pdblines[i]
            final_sections.append(line)
            i += 1

        # Chains selection :
        if chains == 'all' :
            chains_id_list = dict_chains.keys()
        else :
            chains_id_list = sorted(chains)
        pdb_id = self.id
        self.id = list()
        self.path = list()

        # Write the different files
        for chain_id in chains_id_list :
            sub_file_id = pdb_id +  '_' + chain_id
            sub_file_name = 'pdb' + sub_file_id + '.ent'
            sub_file_path = output_dir + sub_file_name
            f = open(sub_file_path, 'w')
            if all_sections :
                f.writelines(initial_sections)
            f.writelines(dict_chains[chain_id])
            if all_sections :
                f.writelines(final_sections)
            f.close()
            self.id.append((pdb_id, chain_id))
            self.path.append(sub_file_path)
        return(dict_chains.keys())


def merge_chains( chains ,  name) :

	os.system('rm %s_scwrfRBD_capped.pdb ' %name )
	for chain in chains :
		os.system("grep 'ACE %s' %s_backboneRBD.pdb  >> %s_scwrfRBD_capped.pdb"  %(chain, name, name ))
		os.system("cat .pdb_%s.ent  >> %s_scwrfRBD_capped.pdb"  %( chain, name ))
		os.system("grep 'NHE %s' %s_backboneRBD.pdb  >> %s_scwrfRBD_capped.pdb"  %(chain, name, name ))


def generate_sequence(sequences) :
    tleap_input='source leaprc.protein.ff14SB \n'
    for seq in range(len(sequences)):
        tleap_input=tleap_input + 'm = sequence { %s } \n ' %(sequences[seq]) + \
        'impose m { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 } { { "N" "CA" "C" "N" -40.0 } { "C" "N" "CA" "C" -60.0 } } \n' + '\n relax m \n savepdb m %s_conf1.pdb  \n'  %( seq )
        #+ 'impose m { 2 3 } { { "C" "CA" "CB" "CG"  60.0 }  }
        #		+ 'impose m { 2 3 } { { "C" "CA" "CB" "CG"  180.0 }  } \nrelax m  \nsavepdb m %s_conf2.pdb\n'  %( seq )\
        #		+ 'impose m { 2 3   } { { "C" "CA" "CB" "CG" 300.0 }  } \nrelax m  \nsavepdb m %s_conf3.pdb \n' %( seq )\
        #		+ 'impose m { 2 3 } {{"C" "CA" "CB" "CG"  120.0 }  } \nrelax m  \nsavepdb m %s_conf4.pdb \n '%( seq )
    tleap_input+='quit \n'
    tleapinput=open('tleap.in','w')
    tleapinput.write(tleap_input)
    tleapinput.close()
    subprocess.check_output('tleap -f tleap.in'  , shell=True)

def mdtraj_aligment(cutoff,input='0_conf1.pdb', ref='out219.pdb'):
    inputtraj=mdtraj.load(input)
    reftraj=mdtraj.load(ref)
    atom_indices=[atom.index for atom in inputtraj.topology.atoms if (atom.name in ['CA' or 'CB' ] and atom.residue.name not in ['ALA', 'NHE', 'ACE'])]
    ref_atom_indices=[atom.index for atom in reftraj.topology.atoms if (atom.name in ['CA' or 'CB' ])]
    inputtraj.superpose(reftraj,atom_indices=None, ref_atom_indices=None )


def mdanalysis_aligment(cutoff,input='0_conf1.pdb', ref='out219.pdb'):
    mobile =mda.Universe(input)
    static=mda.Universe(ref)
    #mobile0 = mobile.select_atoms('name CA or name CB and not resname ALA and not resname ACE and not resname NHE' ).positions - mobile.atoms.center_of_mass()
    #ref0 = ref.select_atoms('name CA or name CB and not resname ALA and not resname ACE and not resname NHE' ).positions - ref.atoms.center_of_mass()
    sele='name CG or name CA or name CB  and not resname ALA and not resname ACE and not resname NHE'

    mobile0 = mobile.select_atoms(sele)
    ref0 = static.select_atoms(sele)
    #print( [ at.resname for at in mobile0.atoms] ,  [ at.resname for at in ref0.atoms])

    align.alignto(mobile, static, select=sele)
    r=rmsd(mobile.select_atoms(sele).positions, static.select_atoms(sele).positions)
    if r < cutoff :
        mobile.atoms.write("%sout.pdb" %input[:-4])
    return r

def pymol_aligment(cutoff,input='out219.pdb'):
	path='../'
	#cmd.fetch('1a38')
	#cmd.load( 'tnf.pdb','tnf' )
	#cmd.load( 'residues_tnfR.pdb' , 'residues_tnfR')
	cmd.load( path +'boxali.pdb','RBD')
	cmd.load( path + input, 'residues_anti')

	rmsd=np.zeros([ len(sequences),3])
	good=[]
	for pep in range(len(sequences)):
	#for pep in range(30,50):
		for conf in ['1'] :#,'2','3','4']
			found=False
			conf_min_val=10
			conf_min=1
			#cmd.sculpt_iterate(pep, state=0, cycles=10)
			#cmd.minimize('1', "MMFF94s", "conjugate gradients", 500,  0.0001, False, 6.0,  8.0)
			cmd.load(path +'%s_conf%s.pdb' %(pep,conf ) )
			ref= 'residues_anti  and (name CA or name CB  )  '

			mobile ='%s_conf%s and (name CA or name CB ) and (  not resname  ALA and not resname  ACE ) ' %(pep,conf)
			print('align ' + mobile  + ' , ', ref)

			blabla =cmd.align( mobile, ref, cycles=20)
			print(blabla)
			#blabla= cmd.super( mobile, ref )
			#cmd.cealign('2lchaind', '1')
			#print (blabla[0])



			if blabla[0] < cutoff  :
					output.write(str('%s :%s '   %('%s_conf%s'  %(pep,conf), blabla[0])))
					#print (process_good('%s_conf%s'  %(pep,conf)) ,  clashes_cutoff  ,'%s_peptide'%(pep) )
					if process_good('%s_conf%s'  %(pep,conf)) < clashes_cutoff:
						conf_min_val = blabla[0]
						conf_min = conf
						found=True

					else:
						cmd.delete('%s_peptide'%(pep))
						cmd.delete('bc_%s_peptide'%(pep))
			if found==True  :
				good.append(str('%s_conf%s'  %(pep,conf_min)))
				file=Flat_File()
				file.read_file('%s_scwrfRBD.pdb' % good[-1])
				chains = file.split_PDBfile_by_chains()
				merge_chains(chains ,  good[-1])
				rmsd[pep, 0]=pep
				rmsd[pep, 1]=conf_min
				rmsd[pep, 2]= conf_min_val

	return good , rmsd



def process_good(i):

		#ABC =  '  ( 1a38 and  chain A ) or ( 1a38 and chain B )'
		ABC='RBD'


		cmd.create( i+ '_2' , i )
		cmd.create( i+ '_3' , i )
		cmd.alter( i ,'chain="D"' )
		cmd.alter( i+'_2' , 'chain="E"' )
		cmd.alter( i+'_3' , 'chain="F"' )

		ref2= 'residues_anti and chain S and  (name CA or name CB or name CG ) '
		ref3= 'residues_anti and chain T and  (name CA or name CB or name CG ) '
		mobile2 ='%s_2 and (name CA or name CB or name CG ) and (not resname  ALA and not resname  ACE  ) ' %(i)
		mobile3 ='%s_3 and (name CA or name CB or name CG ) and (not  resname  ALA and not  resname  ACE  ) ' %(i)

		cmd.align( mobile2, ref2)
		cmd.align( mobile3, ref3)
		#print( ' %s  or  %s_2 or %s ' %(i,i,ABC) )
		#cmd.create  ( '%s_peptide'  %(i),   ' %s  or  %s_2 or %s ' %(i,i,ABC) )
		cmd.create  ( '%s_peptide'  %(i),   ' %s  or  %s_2 or %s_3 or %s ' %(i,i,i, ABC) )
		cmd.create  ( '%s_backbone'  %(i),   ' ((%s  or  %s_2 or %s_3  ) and backbone )' %(i,i,i) )  # or resname  in (ACE NHE)
		cmd.create  ( '%s_backboneRBD'  %(i),   '%s_backbone  or %s' %(i,ABC) )
		cmd.create('bc_%s_peptide'%(i),'%s_peptide'%(i) )


		cmd.set('state', 1)  # only
		cmd.set('sculpt_field_mask', 0x200)   # cSculptaVOID
		#cmd.select('to_be_modelled',  )
		cmd.sculpt_activate('bc_%s_peptide'%(i), '1')
		#strain = cmd.sculpt_iterate('bc_%s_peptide'  %(i), cycles=10)
		#cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
		strain = cmd.sculpt_iterate('bc_%s_peptide'  %(i), cycles=50)
		if strain < clashes_cutoff :
			output.write(str('VDW Strain in structure  %s_peptide: %f\n'   %(i, strain)))
		#show_bumps(, 'bump_check')

			cmd.save  ('%s_peptide.pdb'%(i)  , ' %s_peptide'  %(i))
		cmd.save  ('%s_backboneRBD.pdb'%(i)  , ' %s_peptide'  %(i))
		os.system('Scwrl4 -i %s_backboneRBD.pdb -h -o %s_scwrfRBD.pdb -g %s_scwrfRBD'  %(i,i,i) )

		cmd.delete ( i+'_2' )
		#cmd.show_as('cgo','bc_%s_peptide' %(i) )

		return strain

def coordistance (a,b) :
    return math.sqrt((a[0]- b[0])**2 + (a[1]- b[1])**2 + (a[2]- b[2])**2  )

def clash_score (structure1,structure2 ) :
    u1 = mda.Universe(structure1).atoms.positions
    u2 = mda.Universe(structure2).atoms.positions
    clash=0
    for atom_pos1 in u1 :
        for atom_pos2 in u2 :

            if coordistance ( atom_pos1, atom_pos2) < 2 :
                #print(clash ,coordistance ( atom_pos1, atom_pos2) < 2)
                atom_pos1
                clash=clash+1
    return clash





def positionsidechain(input='1_conf1out.pdb', ref='sameantibodychain169.pdb'):

    static=mda.Universe(ref)

    #mobile0 = mobile.select_atoms('name CA or name CB and not resname ALA and not resname ACE and not resname NHE' ).positions - mobile.atoms.center_of_mass()
    #ref0 = ref.select_atoms('name CA or name CB and not resname ALA and not resname ACE and not resname NHE' ).positions - ref.atoms.center_of_mass()
    dihedrals=[]
    u = mda.Universe(input)
    u.atoms.guess_bonds()
    for dih in u.dihedrals :
        dihe= [ at.name for at in dih.atoms ]
        dihe.sort()
        if dihe == ['C' , 'CA' , 'CB' ,'CG']:
            dihedrals.append(dih)

    for j, dih in enumerate(dihedrals) :
        bvec = dih[2].position -dih[1].position
        ag = dih[1].position
        head = u.select_atoms('resid %s and not backbone and not name H and not name HA '  %dih[0].resid  )
        sele='name CG '
        refs=static.select_atoms(sele).positions
        minr=100

        for i in range(36):
            head.rotateby(i ,  bvec, point=ag)
            r=rmsd(head.select_atoms(sele).positions[0], refs[j])
            r=coordistance (head.select_atoms(sele).positions[0], static.select_atoms(sele).positions[j])
            #sele='not resname ALA and not resname ACE and not resname NHE and not name H*'
            #r=rmsd(u.select_atoms(sele).positions,static.select_atoms(sele).positions)
            if r < minr :
                #u.atoms.write('firstout%s%s.pdb' %( i, j) )
                minr = r
                minangle = i
                new=u.copy()
        u=new
        #u.atoms.write('firstout%s.pdb' %( j) )
    u.atoms.write('firstout.pdb' )
    dihedrals=[]
    for dih in u.dihedrals :
        dihe= [ at.name for at in dih[1:3].atoms ]
        dihe.sort()
        if  ['CB' ,'CG' ]  ==  dihe and 'H' not in [  at.type for at in [dih[0], dih[3] ] ] and 'CA' in[  at.name for at in [dih[0], dih[3] ] ] :
            dihedrals.append(dih)
    for j, dih in enumerate(dihedrals) :
        bvec = dih[2].position -dih[1].position
        ag = dih[1].position
        head = u.select_atoms('resid %s and not backbone and not name H and not name HA  and not name CB and not name HB* '  %dih[0].resid  )
        for i in range(36):
            head.rotateby(i ,  bvec, point=ag)
            #r=rmsd(head.select_atoms(sele).positions[0], refs[j])
            sele='not resname ALA and not resname ACE and not resname NHE and not name H*'
            r=rmsd(u.select_atoms(sele).positions,static.select_atoms(sele).positions)
            if r < minr :
                minr = r
                minangle = i
                new=u.copy()
        head.rotateby(minangle ,  bvec, point=ag)
    u=new
    with mda.Writer('%s-sidechainoriented.pdb' %input[:-4], n_atoms=len(u.atoms)) as w:
        w.write(u.atoms)





sequences=[]
space = [2,3,4]

input=sys.argv[1]

oups=[]
seqs=[]
#seqs.append('TRP LEU ASP LEU'.split())
#seqs.append('LEU ASP LEU'.split())

traj = mdtraj.load(input)
#seqs.append('LEU GLN ARG LEU'.split())
seqs.append([res.name for res  in  traj.topology.residues ])


new_added=True

#mut={ 'ASN':['GLN'] , 'ASP':['GLU'], 'GLN':['ASN'], 'GLU':['ASP'], 'HIS':['TRP'], 'ILE' : ['LEU','VAL'] , 'LEU' :['ILE', 'VAL'], 'PHE' :['TYR'],'SER':['THR'],'THR':['SER'],'TRP' :['HIS'], 'TYR' :['PHE'], 'VAL' :['LEU','ILE']}

#mut={ 'ASN':['GLN'] , 'ASP':['GLU'], 'GLN':['ASN'], 'GLU':['ASP'], 'HIS':['TRP'], 'ILE' : ['LEU','VAL'] , 'LEU' :['ILE', 'VAL'],'SER':['THR'],'THR':['SER'], 'VAL' :['LEU','ILE']}
#mut={ 'ASN’:[‘GLN’] }

'''
while new_added:
	print(seqs)
	new_added=False
	for seq in seqs:
		print( seq)
		for res in range(len(seq)) :
			print( res)
			if seq[res] in mut.keys() :
				for value in mut[seq[res]]:
				 new=copy.deepcopy(seq)
				 new[res]=value
				 if new not in seqs :
				  seqs.append(new)
				  new_added=True
'''

for seq in seqs:
    for j in space :
        for i in space :
            if len(seq) == 4:
                print('lenght sequence = 4 ')
                for k in space :
                   newseq= 'ACE ' +  seq[0]  + ' ' + 'ALA ' * j + ' ' + seq[1] + ' ' +'ALA ' * i +' ' + seq[2] +' ' +  'ALA ' * k + ' ' +seq[3] + ' '   +' NHE'
                   #print(newseq)
                   if newseq not in oups : oups.append(newseq )
            if len(seq) == 3:
               newseq= 'ACE  ALA  ' +  seq[0]  + ' ' + 'ALA ' * j + ' ' + seq[1] + ' ' +'ALA ' * i +' ' + seq[2] + ' ' +  'ALA NHE'
               #print( newseq)
               if newseq not in oups : oups.append(newseq )
for value in oups :
	oups.index(value)
	if len(value.split()) < 14 :
		sequences.append(value)




output=open('out.log', 'w')


generate_sequence(sequences)
'''
for cut in rmsd_cutoff :
    #os.mkdir('cutoff%s'  %cut)
    os.chdir('cutoff%s'  %cut)
    good, rmsd = pymol_aligment(cut)
    os.chdir('../')
'''
good=[]
dirs = os.listdir('.')

cutoff=1.9
for seq in [ file for file in dirs if 'conf1.pdb' in  file ] :
    r = mdanalysis_aligment(cutoff=cutoff, input=seq, ref=input)
    if r < cutoff :
        name=seq[:-4]+'out'
        good.append(name)

        positionsidechain(input=name+'.pdb', ref=input)
print(good)


f= open(input)
header=f.readline()
print(header)

RBDs = [str.split('-')[0] for str in header.split()[1:]]
RBDs = list(dict.fromkeys(RBDs))
print(RBDs )
clashmax=400
for k, struct in enumerate(good) :
    struct=struct+'-sidechainoriented'
    for i,RBD in enumerate(RBDs) :
        if isfile('../RBD/%s.pdb'%RBD) :
            pass
        elif isfile('../RBD/%s.pdb.side'%RBD) :
            print('RBD %s wasnt used initial structure not suitable' %RBD )
            RBDs.remove(RBD)
        else : sys.exit('RBD %s not found !' %RBD)
    for RBD in RBDs :

        print('../RBD/%s.pdb'%RBD , '%s'%struct )
        clash= clash_score('../RBD/%s.pdb'%RBD , '%s.pdb'%struct )
        print('clash score : %s ' %clash)

        if clash < clashmax:
            print('%s with %s passed the class test !' %( RBD  ,struct))
            clashmax=clash

            os.system('grep ATOM ../RBD/%s.pdb  > temp.pdb ' %(RBD) )
            f=open( 'temp.pdb' ,'r' )
            chain = f.readline().split()

            while chain[0] != 'ATOM' :
                chain = f.readline().split()
            print(chain)
            chain=chain[4]
            f.close()
            f=open('%s.pdb' % struct )


            chainpep = f.readline().split()
            while chainpep[0] != 'ATOM' : chainpep = f.readline().split()
            chainpep=chainpep[4]
            f.close()
            '''
            os.system('grep ATOM ../RBD/%s.pdb  > output.pdb ' %(RBD) )
            os.system('grep ATOM %s.pdb  >> output.pdb ' %( struct ) )
            os.system('grep ATOM %s.pdb > outputpep.pdb' %(struct ) )
            os.system('pdb_seg output.pdb > out.pdb')
            out= mda.Universe('out.pdb')
            out.select_atoms('all').write('output.pdb')
            '''
            #print('grep ATOM  %s.pdb | pdb_rplchain -%s:P >> output.pdb ' %( struct, chainpep ))
            os.system('grep HEADER %s |   tr -d \'\n\'  > output%s.pdb ; echo \' %s\'  >> output%s.pdb ' %(input, k, RBD, k) )
            os.system('grep ATOM ../RBD/%s.pdb  |  pdb_rplchain -%s:A  >> output%s.pdb ' %(RBD, chain ,k) )
            os.system('grep ATOM  %s.pdb | pdb_rplchain -%s:P >> output%s.pdb ' %( struct , chainpep ,k) )
            os.system('grep ATOM  %s.pdb | pdb_rplchain -%s:P > outputpep%s.pdb' %( struct, chainpep,k ) )
            final= struct

'''
def loadstructures(mypath='../'):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[-4:]=='.pdb' and '_ali' in f]
    atoms={}
    for f in onlyfiles :
        print(f)

        atoms.update({f[:-4]:mda.Universe(mypath+f)})
    return atoms

final= mda.Universe(struct +'.pdb')
mutationspoints=final.select_atoms('resname ALA and name CB').positions
resindex = final.select_atoms('resname ALA and name CB').resids
print(resindex)
sidechains=  loadstructures()
print()
for j, pos in enumerate(mutationspoints) :
    print(pos)
    for i , (sidechain, traj) in  enumerate(sidechains.items()) :
        new_pos= traj.select_atoms('name CB').positions
        if len(new_pos)>0 :
            new_pos=new_pos[0]
            #print(coordistance ( pos, new_pos))
            if coordistance ( pos, new_pos) <3 :
                print('mutation possible'  + str(resindex[j]) +  traj.select_atoms('name CB').resnames[0])
    sidechains
'''






#cmd.group('initial' , '*_conf*')
#cmd.group('to_check' , '*_conf*peptide')
'''
for sequence in good  :
	file=Flat_File()
	file.read_file('%s_scwrfRBD.pdb' % sequence)
	chains = file.split_PDBfile_by_chains()
	merge_chains(chains , sequence)
rmsd.sort(axis=0)
print ( str(rmsd))

output.close()
cmd.disable()
'''



# for visualising our result
