from Bio.PDB import *
import os
import argparse
import mdtraj as md
import os
import subprocess
import time
import shutil

dir='../'
tleapheader='set default PBradii mbondi2\n\
source leaprc.protein.ff14SB \n\
source leaprc.water.tip4pew \n\
loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/S51.lib \n\
loadamberparams /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/parmchk2_all.frcmod \n\
loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/S52.lib \n\
loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/R52.lib \n\
loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/R51.lib \n'




def remove_res(inputfile, outputfile,  delchains=['D', 'E' , 'F' ],  posid=0):
    "remove amino acid at position posid : 0 first, -1 last"

    chains=split_PDBfile_by_chains( inputfile)

    for chain in delchains :
        system = md.load_pdb('.pdb_%s.ent' %(chain))
        topology =system.topology
        peptide = topology.chain(0)
        firstres=peptide.residue(0).resSeq
        lastres=peptide.residue(-1).resSeq
        if posid ==0 :
            firstres=peptide.residue(1).resSeq
            remove_resname_pdb( '.pdb_%s.ent' %(chain) , 'del.pdb_%s.ent' %chain, posid , firstres)

        if posid ==-1 :
            lastres=peptide.residue(-2).resSeq
            remove_resname_pdb( '.pdb_%s.ent' %(chain) , 'del.pdb_%s.ent' %chain , posid ,lastres )
        fill_tleap('del.pdb_%s.ent' %(chain)   , 'delf.pdb_%s.ent' %(chain)  , firstres,chain )
        pdbs2merge=[]
    for i in  chains:
        if i in delchains :
            pdbs2merge.append(  'delf.pdb_%s.ent' %i  )
        else : pdbs2merge.append( '.pdb_%s.ent' %i)
        mergepdbs(pdbs2merge,outpoutname= outputfile )


def remove_resname_pdb(input, output, posid, resid) :
    "Mdtraj just renumber , pdb4amber sometimes crash for no reason... just came up with my oow pdb parser !\
    takes name of the input, name of the output (SINGLE CHAIN PDB ! and mutations (the fonction will apply mutations on the chain) )"
    origin = open( input, 'r')
    originlines=origin.readlines()
    mutated = open( output, 'w')
    for line in originlines :
        linelist= list (line)
        if posid ==0 and 'ACE' in line :     pass

        elif  posid ==0 and line.split()[0] == 'ATOM' and   int(line[23:26]) == resid :
            if (line[13:15] in  ['C ', 'O '  ]  )  :
                linelist[17:20]=list ('ACE')
                line = "".join(linelist)
                mutated.writelines(line)

        elif posid ==-1 and 'NHE' in line : pass

        elif  posid ==-1 and line.split()[0] == 'ATOM' and   int(line[23:26]) == resid:
            if (line[13:15] in  ['N ', 'H '  ]  )  :
                linelist[17:20]=list ('NHE')
                line = "".join(linelist)
                mutated.writelines(line)


        else : mutated.writelines(line)
    origin.close()
    mutated.close()




def fill_tleap(  inputfile,outputfile,firstres, chain) :
            os.system('grep -v H.3  %s | grep -v H.1 | grep -v CD\\ \\ ILE | grep -v HD2\\ ILE > no3H-no1H-%s' %(  inputfile,   inputfile))
            tleapin = open ('tleap_chain.in', 'w' )
            #tleapin.writelines(tleapheader)
            #stapledpeptides.save_pdb('828_conf1_scwrftnf_capped_mod%s.pdb' %i)

            tleapin.writelines('source leaprc.protein.ff14SB \n' )
            tleapin.writelines('a = loadpdb no3H-no1H-%s\n' %inputfile )
            tleapin.writelines('select a\n')
            tleapin.writelines('relax a\n')
            tleapin.writelines('savepdb a out.pdb  \n')

            #os.system('grep ATOM %s_stapled%s-tnf.pdb > temp ; mv  temp %s_stapled%s-tnf.pdb ' %(file[:-4],i,file[:-4],i))
            tleapin.writelines('quit\n')
            tleapin.close()
            subprocess.check_output('tleap -f tleap_chain.in' ,shell=True)
            #leave som time to tleap to close all files properly
            time.sleep(3)


            os.system( 'pdb_reres -%s  out.pdb> out-res.pdb' % (firstres))
            p=open('out-res.pdb' ,'r')
            to_rechain= p.readlines()
            p.close()
            p=open('out-res-chain.pdb','w')
            for l in to_rechain:
                if 'ATOM' in  l : p.writelines('%s%s%s' %(l[:21] , chain,l[22:] ))
                else :  p.writelines(l)
            p.close()
            os.system( 'pdb_reatom -1 out-res-chain.pdb >  %s'%(outputfile))

def add_staple(inputfile, chain='D',  stapledres1='S51',stapledres2='S52', res2ditch=['ALA', 'SER'] , allbut=False , posid=0):
    "if posid=0 put a staple everywhere thare is a res at position i and i+4 else put a staple at posid "

    chains=split_PDBfile_by_chains( inputfile)
    system = md.load_pdb('.pdb_%s.ent' %(chain))
    topology =system.topology
    peptide = topology.chain(0)
    firstres=peptide.residue(0).resSeq
    if posid==0:
        for i  in range(peptide.n_residues-4)  :

            if (  res2ditch == 'all' or ( peptide.residue(i).name in res2ditch  and  peptide.residue(i+4).name in res2ditch ) ) and (allbut == False or (peptide.residue(i).resSeq not  in allbut and    peptide.residue(i+4).resSeq not  in allbut )) :
                print( 'staple' , i )

                change_resname_pdb( '.pdb_%s.ent' %(chain),  'tobefilled.pdb', [[chain , peptide.residue(i).resSeq,  stapledres1 ] ,[chain ,  peptide.residue(i+4).resSeq, stapledres2] ])
                #fill_staple_tleap( 'tobefilled.pdb',  '%s-stp%s.pdb' %  (inputfile[:-9], peptide.residue(i).resSeq) , stapledres1,stapledres2 )
                fill_staple_tleap( 'tobefilled.pdb',  '%s-stp%s.pdb'  % (inputfile[:-9], peptide.residue(i).resSeq) , stapledres1,stapledres2 ,firstres, chain )





                pdbs2merge=[]
                for j in  chains:
                    if j == chain :
                        pdbs2merge.append(  '%s-stp%s.pdb' % (inputfile[:-9], peptide.residue(i).resSeq))
                    else : pdbs2merge.append( '.pdb_%s.ent' %j)
                    mergepdbs(pdbs2merge,outpoutname= '%s-stp%s.pdb' % (inputfile[:-4],peptide.residue(i).resSeq ))
    else :
            i = posid
            change_resname_pdb(  inputfile , 'tobefilled.pdb', [[chain , i , stapledres1 ] ,[chain , i+4 , stapledres2  ] ])
            fill_staple_tleap( 'tobefilled.pdb',  '%s-stp%s.pdb' % (inputfile[:-4],i) ,stapledres1,   stapledres2)
            pdbs2merge=[]
            for i in  chains:
                if i == chain :
                    pdbs2merge.append(  '%s-stp%s.pdb' % (inputfile[:-4],i))
                else : pdbs2merge.append( '.pdb_%s.ent' %i)
                mergepdbs(pdbs2merge,outpoutname= '%s-stp%s.pdb' % (inputfile[:-4], i ))

def add_fluo(inputfile, chain='D',  res='DYE', res2ditch=['ALA', 'SER'] , allbut=False , posid=0):
    "if posid=0 put a staple everywhere there is a res at position i and i+4 else put a staple at posid "

    chains=split_PDBfile_by_chains( inputfile)
    system = md.load_pdb('.pdb_%s.ent' %(chain))
    topology =system.topology
    peptide = topology.chain(0)
    firstres=peptide.residue(0).resSeq
    if posid==0:
        for i  in range(peptide.n_residues-4)  :
            print(peptide.residue(i).name, peptide.residue(i+4).name, peptide.residue(i).resSeq)

            if (  res2ditch == 'all' or ( peptide.residue(i).name in res2ditch  )) and (allbut == False or (peptide.residue(i).resSeq not  in allbut )) :

                change_resname_pdb( '.pdb_%s.ent' %(chain),  'tobefilled.pdb', [[chain , peptide.residue(i).resSeq,  stapledres1 ] ,[chain ,  peptide.residue(i+4).resSeq, stapledres2] ])
                #fill_staple_tleap( 'tobefilled.pdb',  '%s-stp%s.pdb' %  (inputfile[:-9], peptide.residue(i).resSeq) , stapledres1,stapledres2 )
                fill_staple_tleap( 'tobefilled.pdb',  '%s-fluo%s.pdb'  % (inputfile[:-9], peptide.residue(i).resSeq) , stapledres1,stapledres2 ,firstres, chain )
                pdbs2merge=[]
                for j in  chains:
                    if j == chain :
                        pdbs2merge.append(  '%s-fluo%s.pdb' % (inputfile[:-9], peptide.residue(i).resSeq))
                    else : pdbs2merge.append( '.pdb_%s.ent' %j)
                mergepdbs(pdbs2merge,outpoutname= '%s-fluo%s-full.pdb' % (inputfile[:-4],peptide.residue(i).resSeq ))
    else :
            i = posid
            change_resname_pdb( '.pdb_%s.ent' %(chain)  , 'tobefilled.pdb', [ [ chain,  i , res ]  ])
            fill_fluo_tleap( 'tobefilled.pdb','%s-fluo%s.pdb' % (inputfile[:-4],i),res, firstres=firstres )
            pdbs2merge=[]
            for j in  chains:
                if j == chain :
                    pdbs2merge.append(  '%s-fluo%s.pdb' % (inputfile[:-4],i))
                else :
                    print( '.pdb_%s.ent' %j)
                    pdbs2merge.append( '.pdb_%s.ent' %j)
            mergepdbs(pdbs2merge,outpoutname= '%s-fluo%s-full.pdb' % (inputfile[:-4] ,i))

def fill_fluo_tleap(  inputfile,outputfile, RES1, firstres=1, chain='D') :
            os.system('grep -v H.3  %s | grep -v H.1 > no3H-no1H-%s' %(  inputfile,   inputfile))
            tleapin = open ('tleap_chain.in', 'w' )
            #tleapin.writelines(tleapheader)
            #stapledpeptides.save_pdb('828_conf1_scwrftnf_capped_mod%s.pdb' %i)

            tleapin.writelines('source leaprc.protein.ff14SB \n' )
            tleapin.writelines('loadoff Dye.lib \n') # %(RES1)
            tleapin.writelines('loadamberparams dye.frcmod \n' )
            tleapin.writelines('a = loadpdb no3H-no1H-%s\n' %inputfile )
            tleapin.writelines('savepdb a out.pdb  \n')
            #os.system('grep ATOM %s_stapled%s-tnf.pdb > temp ; mv  temp %s_stapled%s-tnf.pdb ' %(file[:-4],i,file[:-4],i))
            tleapin.writelines('quit\n')
            tleapin.close()
            os.system('tleap -f tleap_chain.in')
            #leave som time to tleap to close all files properly
            time.sleep(3)
            os.system( 'grep -v TER  out.pdb> out-noter.pdb')

            os.system( 'pdb_reres -%s  out-noter.pdb> out-res.pdb' % (firstres))
            p=open('out-res.pdb' ,'r')
            to_rechain= p.readlines()
            p.close()
            p=open('out-res-chain.pdb','w')
            for l in to_rechain:
                if 'ATOM' in  l : p.writelines('%s%s%s' %(l[:21] , chain,l[22:] ))
                else :  p.writelines(l)
            p.close()
            os.system( 'pdb_reatom -1 out-res-chain.pdb >  %s'%(outputfile))

            #cmd='file=%s_stapled%s-tnf.pdb ; grep -v HH31\\ ACE $file > temp ;  grep -v HH32\\ ACE temp > $file ;  grep -v HH33\\ ACE $file > temp ; mv temp $file  '%( file[:-4],i)
            #print(cmd)
            #os.system(cmd)
            #cmd='file=%s_stapled_%i.pdb; do  grep -v H1\\ \\ ACE $file > temp ;  grep -v H2\\ \\ ACE temp > $file ;  grep -v H3\\ \\ ACE $file > temp ; mv temp $file ;done  '%(inputfile[:-4],i)
            #os.system(cmd)





def fill_staple_tleap(  inputfile,outputfile, RES1,RES2 ,firstres=1, chain='D') :
            os.system('grep -v H.3  %s | grep -v H.1 > no3H-no1H-%s' %(  inputfile,   inputfile))
            tleapin = open ('tleap_chain.in', 'w' )
            #tleapin.writelines(tleapheader)
            #stapledpeptides.save_pdb('828_conf1_scwrftnf_capped_mod%s.pdb' %i)

            tleapin.writelines('source leaprc.protein.ff14SB \n' )
            tleapin.writelines('loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/%s.lib \n' %(RES1) )
            tleapin.writelines('loadamberparams /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/parmchk2_all.frcmod \n' )
            tleapin.writelines('loadoff /home/marie/Desktop/newerparam-again/SS-bisubstitued-stapled-cis/%s.lib \n' %(RES2) )

            tleapin.writelines('a = loadpdb no3H-no1H-%s\n' %inputfile )
            tleapin.writelines('bond a.%s.CY a.%s.CY \n' %(RES1,RES2 ))
            tleapin.writelines('select a\n')
            tleapin.writelines('relax a\n')
            tleapin.writelines('savepdb a out.pdb  \n')

            #os.system('grep ATOM %s_stapled%s-tnf.pdb > temp ; mv  temp %s_stapled%s-tnf.pdb ' %(file[:-4],i,file[:-4],i))
            tleapin.writelines('quit\n')
            tleapin.close()
            os.system('tleap -f tleap_chain.in')
            #leave som time to tleap to close all files properly
            time.sleep(3)


            os.system( 'pdb_reres -%s  out.pdb> out-res.pdb' % (firstres))
            p=open('out-res.pdb' ,'r')
            to_rechain= p.readlines()
            p.close()
            p=open('out-res-chain.pdb','w')
            for l in to_rechain:
                if 'ATOM' in  l : p.writelines('%s%s%s' %(l[:21] , chain,l[22:] ))
                else :  p.writelines(l)
            p.close()
            os.system( 'pdb_reatom -1 out-res-chain.pdb >  %s'%(outputfile))
            #cmd='file=%s_stapled%s-tnf.pdb ; grep -v HH31\\ ACE $file > temp ;  grep -v HH32\\ ACE temp > $file ;  grep -v HH33\\ ACE $file > temp ; mv temp $file  '%( file[:-4],i)
            #print(cmd)
            #os.system(cmd)
            #cmd='file=%s_stapled_%i.pdb; do  grep -v H1\\ \\ ACE $file > temp ;  grep -v H2\\ \\ ACE temp > $file ;  grep -v H3\\ \\ ACE $file > temp ; mv temp $file ;done  '%(inputfile[:-4],i)
            #os.system(cmd)


def prep_filexfold(inputfile,chains2refine ):
    chains=split_PDBfile_by_chains( inputfile )
    pdbs2merge  = []
    for i in  chains:
            if i in  chains2refine :
                pdbs2merge.append( '.pdb_%s.ent' %(str(i)) )
            else :
                shutil.copy('.pdb_%s.ent' %(str(i)),'copy.pdb_%s.ent' %(str(i)) )
                pdbs2merge.append('copy.pdb_%s.ent' %i )
    os.system("grep CYS  output0.pdb | awk  '{print $5 $6}'  > cys.txt " )
    f=open('cys.txt')
    cysres = 'C' + ',C'.join(list(dict.fromkeys(f.readlines())))
    f.close()
    f=open('cys.txt' ,'w')
    f.writelines(cysres.replace('\n','')+'\n')
    f.close()
    os.system('~/Utilities/foldx5Linux64/foldx --command=RepairPDB --fix-residues-file=cys.txt  --pdb=' +inputfile)
    os.system("sed -i '/^$/d' %s_Repair.pdb"   %inputfile[:-4] )
    chains=split_PDBfile_by_chains( inputfile[:-4] + '_Repair.pdb' )
    mergepdbs(pdbs2merge ,outpoutname=inputfile[:-4] + '_prep.pdb')


def split_PDBfile_by_chains( path , id='', output_dir = '.', chains = 'all', all_sections = True ) :
    ''' Split a pdb file in different pdb files by chains. data is a list of
    pdb file lines. chains must be a list of PDB ids (e.g. ['A', 'B'])
    '''
    if path != '.' :
        os.system('grep -v SSBOND %s >temp.pdb' %  path)
        f = open('temp.pdb', "r")

        pdblines = f.readlines()
        f.close()

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
        possible_sections = ['ATOM  ', 'ANISOU', 'HETATM']
        if line[0:6]in possible_sections:
            chain_id = line[21]

            if not(chain_id in dict_chains) :
                dict_chains[chain_id] = [line]
            else :
                dict_chains[chain_id].append(line)
            i += 1
        elif 'TER' in line :i += 1
        elif 'SSBOND' in line :i += 1
        else :
            break
    while i < len(pdblines) :
        line = pdblines[i]
        final_sections.append(line)
        i += 1

    # Chains selection :
    if chains == 'all' :
        chains_id_list = dict_chains.keys()

    else :
        chains_id_list = sorted(chains)
    pdb_id = id
    id = list()
    path = list()

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
        id.append((pdb_id, chain_id))
        path.append(sub_file_path)
    return(dict_chains.keys())




def mergepdbs(pdbs,outpoutname='temp'):
    # bash call faster than python :-(
    os.system('rm %s' %(outpoutname))
    for pdb in pdbs :
            os.system('grep HEADER   %s  >> %s' %(pdb,outpoutname))

    for pdb in pdbs :
            os.system('grep -v TER   %s | grep -v END |  grep -v HEADER   >> %s' %(pdb,outpoutname))
            #os.system('cat %s >> %s' %(pdb,outpoutname))
            os.system('echo TER >> %s' %(outpoutname))

def change_resname_pdb(input, output, mutations) :
    "Mdtraj just renumber , pdb4amber sometimes crash for no reason... just came up with my oow pdb parser !\
    takes name of the input, name of the output (SINGLE CHAIN PDB ! and mutations (the fonction will apply mutations on the chain) )"
    origin = open( input, 'r')
    originlines=origin.readlines()
    mutated = open( output, 'w')
    for line in originlines :
        #print( int(line[23:26]) )
        if  line.split()[0] == 'ATOM' and   int(line[23:26])  in [int(mut[1]) for mut in mutations] :
            mut= mutations[[int(mut[1]) for mut in mutations].index(  int(line[23:26])) ]
            #print('MUTATION' ,mut )
            linelist= list (line)
            if (line[13:15] in  ['C ','N ','O ','CA','H ' , 'HA' , 'CB'  ] and mut[2]!='DYE' ) or line[13:15] in  ['C ','N ','O ','CA'] :

                linelist[17:20]=list (mut[2])
                line = "".join(linelist)
                mutated.writelines(line)
        else : mutated.writelines(line)
    origin.close()
    mutated.close()



def mutate_residues( filename, mutations, stemout='out') :
    print('stemout begin')
    print(stemout)
    chains = split_PDBfile_by_chains(filename)
    print(filename)
    print( chains)
    if mutations :
        mutchains=[]


        for mut in mutations :

            if mut[0] not in mutchains:
                mutchains.append( mut[0])


        # list of list chain residus_nummer and new_residus
        for j in mutchains:
            change_resname_pdb('.pdb_%s.ent' %(j), 'mut.pdb_%s.ent' %(str(j)), mutations)

        pdbs2merge=[]
        print( chains)
        for i in  chains:
            if i in mutchains :
                pdbs2merge.append( 'mut.pdb_%s.ent' %i)
            else :
                pdbs2merge.append( '.pdb_%s.ent' %i)
                #system = md.load_pdb('.pdb_%s.ent' %(i))
                #firstres=system.topology.chain(0).residue(0).resSeq
                shutil.copy('.pdb_%s.ent' %(i), 'copy.pdb_%s.ent' %(i))

                #os.system('~/Utilities/Maestro_2017-3_Linux-x86_64_Academic/utilities/prepwizard  original.pdb_%s.ent prep.pdb_%s.ent -noepik   -disulfides   -noimpref' %(i,i))
        print('PDB2merge 1:')
        print( pdbs2merge)
        mergepdbs(pdbs2merge,outpoutname='temp2.pdb')
        ###  Scwrl4 remove all unatural amino acids !!!!! need to reconstruct pdb
        subprocess.check_output('Scwrl4  -i  temp2.pdb  -o tempScwrl4.pdb',shell=True)
        split_PDBfile_by_chains( 'tempScwrl4.pdb')
        # Now mut.x contains not full side_chains & .x contain no unatural a.a. but full sidecahins Let's make a file that isn't a mess
        for  chain  in mutchains:
            origin = open('mut.pdb_%s.ent' %(chain) , 'r')
            originlines=origin.readlines()

            originformutres = open('.pdb_%s.ent' %(chain) , 'r')
            originlinesformutres=originformutres.readlines()

            output = open('%s_%s.pdb' %(stemout,str(chain)) , 'w')
            fullmutatedres=[]


            for line in originlines :
                #print( int(line[23:26]) )
                if line.split()[0] == 'ATOM'  and  int(line[23:26])  in [int(mut[1]) for mut in mutations if mut[0]==chain ] :
                    if line[13:15] == 'CA':  #do it once only
                        for linemutres in originlinesformutres :
                            if linemutres.split()[0] == 'ATOM' and len(linemutres) > 5 and   int(linemutres[23:26]) ==  int(line[23:26]) :
                                #print(linemutres)
                                output.writelines(linemutres)


                else :output.writelines(line)
            origin.close()
            output.close()
            originformutres.close()
        pdbs2merge=[]
        for i in  chains:
            if i in mutchains :
                pdbs2merge.append( '%s_%s.pdb' %(stemout,str(i)) )
            else :
                pdbs2merge.append('copy.pdb_%s.ent' %i )
        print('PDB2merge:')
        print( pdbs2merge)
        print(chains)
        print(stemout)
        mergepdbs(pdbs2merge,outpoutname='%s_full.pdb' %(stemout))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='mutate residus and produce pdbs files ready for alchemical calculations')
    parser.add_argument('-m',  type=str, default=False , help='mutations in a txt file  each line must include 2 columns : chain  resnumber new_residue_type')
    parser.add_argument('-f',  type=str,    help='initialpdb file')
    parser.add_argument('-s',  type=str,    help='produce stapled peptide on demand enter all for all position or a residue id')

    args = parser.parse_args()
    if args.m :
        mut= open (args.m)
        mutations=[]
        for line in mut.readlines() :
            mutations.append(line.split())
        mut.close()

        mutate_residues( args.f,mutations)
    else : print(' no mutation file given (-m argument), not doing any mutations )')

    if args.s :
        if args.s == 'all' :
            posid=0
        else:  posid=args.s
        add_staple( args.f, chain='D', stapledres1='S51',stapledres2='S52' ,posid=posid)
