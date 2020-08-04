#!/usr/bin/env python3

import os
import math
import argparse


def file_size_check(path,fsize=10):
    """Check file's existence and size (MB)"""

    log = {'nice':True,}
    try:
        sizetmp = os.stat(path).st_size
        if sizetmp/1024/1024 > fsize:
            log['nice'] = False
            log['info'] = 'Error: the file size is far larger than %f MB' % fsize,
            
    except IOError:
        log['nice'] = False
        log['info'] = 'Error : cannot open the file!\n' + \
                      'Error : ' + path
    return log



def file_gen_new(fname,fextend='txt',foriginal=True):
    """Return unique new file name in the current path to avoid overriding old files"""
     
    if foriginal is True:
        try:
            f = open(fname + '.' + fextend)
            f.close()
        except IOError:
            return fname + '.' + fextend
     
    i = 1
    filename = fname
    while True:
        fname = filename + '_' + str(i) + '.' + fextend
        try:
            f = open(fname)
            f.close()
        except IOError:
            break
        i += 1
   
    return fname



class Overlap(object):
    def __init__(self,*args,**kwargs):
        self.periodic_table = { 3:'Li', 4:'Be', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 19:'K',
                                20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe',
                                27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As',
                                34:'Se', 37:'Rb', 38:'Sr', 39:'Y', 40:'Zr', 41:'Nb', 42:'Mo',
                                44:'Ru', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb',
                                52:'Te', 55:'Cs', 56:'Ba', 72:'Hf', 73:'Ta', 74:'W', 75:'Re',
                                76:'Os', 78:'Pt', 79:'Au', 82:'Pb', 83:'Bi', 84:'Po', 85:'At',
                                1:'H', 5:'B', 17:'Cl', 7:'N', 8:'O', 9:'F', 16:'S', 6:'C',
                                35:'Br', 53:'I', 15:'P' }
        
        if 'file' in kwargs and kwargs['file'] is not None:
            self.file = kwargs['file']
            log = file_size_check(self.file,500)
            if not log['nice']:
                print(log['info'])
                exit()
        else:
            print('Error: the parameter file is missing')
            exit()


        if 'fname' in kwargs and kwargs['fname'] is not None:
            self.fname = kwargs['fname']
        else:
            self.fname = 'Overlap'


        if 'pbc' in kwargs and kwargs['pbc'] is False:
            self.pbc = False
        else:
            self.pbc = True


        if 'fixatoms' in kwargs and kwargs['fixatoms'] is True:
            self.fixatoms = True
        else:
            self.fixatoms = False


        if 'distance' in kwargs and kwargs['distance'] is not None:
            try:
                self.distance = float(kwargs['distance'])
            except ValueError:
                print('Error: the parameter distance has to be a number')
                exit()
        else:
            self.distance = 5


        if 'finput_format' in kwargs and kwargs['finput_format'] is not None:
            self.finput_format = kwargs['finput_format'].lower()
            if self.finput_format not in ['pdb','gro','xyz','txt']:
                print('Error: currently the Forced input file format is not supported')
                print(kwargs['finput_format'])
                exit()
        else:
            self.finput_format = None


        if 'fout_format' in kwargs and kwargs['fout_format'] is not None:
            self.fout_format = kwargs['fout_format'].lower()
            if self.fout_format not in ['pdb','gro','xyz']:
                print('Error: the out file format is not supported')
                print(kwargs['fout_format'])
                exit()
        else:
            self.fout_format = 'xyz'

        if self.finput_format is None:
            self.finput_format = self.file.split('.')[-1].lower()
            if self.finput_format not in ['pdb','gro','xyz','txt']:
                print('Error: the input file format is not supported')
                print(self.file)
                exit()

        if self.finput_format == 'txt': self.pbc = False


        # Processing the input file
        # self.profile ( residueChainID, residueName, atomType, cor_x, cor_y, cor_z )
        self.profile = []
        if self.finput_format == 'pdb':
            self.profile_pdb()
        elif self.finput_format == 'gro':
            self.profile_gro()
        elif self.finput_format == 'xyz':
            # Take care of Special Case: xyz file
            self.profile_xyz()
        elif self.finput_format == 'txt':
            # Take care of Special Case: txt file
            self.profile_txt()
        else:
            pass
        
        if len(self.profile) == 0:
            print('Error: no inputs were found')
            exit()

        # Analyzing
        #
        # Final Processed Attributes Format
        #
        # self.corlist: 4D
        #           [     [     [    [ x, y, z ]     ]    ]    ]     
        #                res   mol  atom  
        #
        # self.atomtype: 2D
        #           [     [ atype1, atype2, atype3, ... ],        ]
        #
        # self.reslist: 1D
        #           [     residue_name_1,   residue_name_2        ]
        self._profile()

        # Check Valid, for same residue, the number of molecule's atoms should be the same
        for res in self.corlist:
            t = [ len(mol) for mol in res ]
            if len(set(t)) != 1:
                print('Error: the input file is broken')
                exit()

        if self.pbc:
            tsum = 0
            for res in self.corlist:
                if len(res) != 0:
                    tsum += 1
            if tsum == 0:
                print('Error: the system is not fully equilibrated ')
                print('   + : all the molecules are in their biggest chain length stretchs')
                print('   + : this may due to the size of system is so small')
                print('   + : which caused all molecules inside it was excluded!')
                print('   + : One suggestion is to turn off the pbc, setting it to False')
                exit()

        self.custom_atomtype()
        if self.fixatoms:
            self.setorigin()
        self.corlist = self.overlap(self.corlist)



    def profile_pdb(self):
        with open(self.file,mode='rt') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break
                else:
                    if len(line) > 6 and ( line[:4] == 'ATOM' or line[:6] == 'HETATM' ):
                        try:
                            ls = []
                            ls.append(int(line[22:26]))
                            ls.append(line[17:20])
                            ls.append(line[12:16])
                            ls.append(float(line[30:38]))
                            ls.append(float(line[38:46]))
                            ls.append(float(line[46:54]))
                            self.profile.append(ls)
                        except ValueError:
                            print('Warning: the input line is wrong')
                            print(line)


    
    def profile_gro(self):
        # Note: unit in nanometer, should multiple by 10
        with open(self.file,mode='rt') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break
                else:
                    if len(line) >= 44:
                        try:
                            ls = []
                            ls.append(int(line[0:5]))
                            ls.append(line[5:10])
                            ls.append(line[10:15])
                            ls.append(float(line[20:28])*10)
                            ls.append(float(line[28:36])*10)
                            ls.append(float(line[36:44])*10)
                            self.profile.append(ls)
                        except ValueError:
                            print('Warning: the input line is wrong')
                            print(line)

    
    def profile_txt(self):
        #Special case, need to add residueChainID and residueName
        with open(self.file,mode='rt') as f:
            file = f.readlines()
            
        i = 0
        count = 1
        prolist = []
        while i < len(file):
            line = file[i]
            if line.lower().find('accept') != -1:
                j = i + 1
                while j < len(file):
                    if len(file[j]) == 0:
                        count += 1
                        break
                    ltmp = file[j].split()
                    if len(ltmp) == 0:
                        count += 1
                        break
                    try:
                        if len(ltmp) != 4:
                            raise ValueError
                        ls = []
                        ls.append(count)
                        ls.append(int(ltmp[0]))
                        ls.append(float(ltmp[1]))
                        ls.append(float(ltmp[2]))
                        ls.append(float(ltmp[3]))
                        prolist.append(ls)
                    except ValueError:
                        print('Warning: the input line is wrong')
                        print(file[j])
                    j += 1                  
                i = j
            else:
                i += 1

        # split to Reference, FEP_1, FEP_2
        i = 0
        reslist = []
        while i < len(prolist):
            j = i
            ls = []
            ndx = prolist[i][0]
            while j < len(prolist):
                if prolist[j][0] == ndx:
                    ls.append(prolist[j][1:])
                else:
                    break
                j += 1
            reslist.append(ls)
            i = j

        if len(reslist) % 3 != 0:
            print('Warning: the input file is broken, truncation will happen')

        for i, mol in enumerate(reslist):
            if i % 3 == 0:
                for cor in mol:
                    self.profile.append( [i, 'REF'] + cor )
            elif i % 3 == 1:
                for cor in mol:
                    self.profile.append( [i, 'EP1'] + cor )
            else:
                for cor in mol:
                    self.profile.append( [i, 'EP2'] + cor )



    def profile_xyz(self):
        pass



    def _profile(self):

        # split all parameters to molecule level
        prolist = []
        count = 0
        pbc_min = min(self.profile[0][3:])
        pbc_max = max(self.profile[0][3:])
        while count < len(self.profile):  
            i = count
            ls = []
            resndx = self.profile[count][0]
            while i < len(self.profile):
                if pbc_min > min(self.profile[i][3:]): pbc_min = min(self.profile[i][3:])
                if pbc_max < max(self.profile[i][3:]): pbc_max = max(self.profile[i][3:])
          
                if self.profile[i][0] == resndx:
                    ls.append(self.profile[i][1:])
                else:               
                    break
                i += 1
            prolist.append(ls)           
            count = i

        # conclude to residue level -- only the molecule's FIRST residueName is used
        syslist = []
        reslist = []
        count = 0
        while count < len(prolist):
            resndx = prolist[count][0][0]
            i = count
            ls = []
            while i < len(prolist):
                if prolist[i][0][0] == resndx:
                    ls.append(prolist[i])
                else:
                    break
                i += 1
            syslist.append(ls)
            reslist.append(resndx)
            count = i

        # combine residues -- if it has any
        if len(reslist) != 1 and len(reslist) != len(set(reslist)):
            # Idea: unique_ndx_list: reflist [0 , 1,  4, ...]
            #       repeat_ndx_list: replist [  2,  3,   ...]
            i = 0
            replist = []
            reflist = []
            while i < len(reslist):
                j = 0
                bo = True
                while j < i:
                    if reslist[j] == reslist[i]:
                        bo = False
                        break
                    j += 1
                if bo:
                    reflist.append(i)
                else:
                    replist.append(i)
                i += 1
            
            unisyslist = []
            proreslist = []
            for resid in reflist:
                res = syslist[resid]
                for ndx in replist:
                    if reslist[resid] == reslist[ndx]:
                        for mol in syslist[ndx]:
                            res.append(mol)
                unisyslist.append(res)
                proreslist.append(reslist[resid])
        else:
            proreslist = reslist
            unisyslist = syslist

        # for residue's each molecule, Remove all the bad PBC molecules
        if self.pbc:
            prosyslist = []
            pbc_half = (pbc_max - pbc_min) / 2
            for res in unisyslist:
                ls = []
                for mol in res:
                    rxlist = []
                    rylist = []
                    rzlist = []
                    for atom in mol:
                        rxlist.append(atom[2])
                        rylist.append(atom[3])
                        rzlist.append(atom[4])
                    if max(rxlist) - min(rxlist) < pbc_half and \
                       max(rylist) - min(rylist) < pbc_half and \
                       max(rzlist) - min(rzlist) < pbc_half:
                        ls.append(mol)
                prosyslist.append(ls)
        else:
            prosyslist = unisyslist


        self.corlist = []
        self.reslist = []
        self.atomtype = []
        for ndx, res in enumerate(prosyslist):
            if len(res) != 0:
                lp = []
                for mol in res:
                    ls = []
                    for atom in mol:
                        ls.append(atom[2:])
                    lp.append(ls)
                self.corlist.append(lp)

                self.reslist.append(proreslist[ndx])

                ls = []
                for atom in res[0]:
                    ls.append(atom[1])
                self.atomtype.append(ls)
            else:
                print('Warning: the residue < {:} > is excluded from processing...'.format(proreslist[ndx]))



    def custom_atomtype(self):

        # Guess real atomType
        if self.finput_format == 'txt':
            keylist = [ i for i in self.periodic_table ]
            for i, mol in enumerate(self.atomtype):
                for j, atom in enumerate(mol):
                    if atom in keylist:
                        self.atomtype[i][j] = self.periodic_table[atom]
            return
        else:
            pt = [ self.periodic_table[i] for i in self.periodic_table ]
            for i, mol in enumerate(self.atomtype):
                for j, atom in enumerate(mol):                    
                    for name in pt:
                        if atom.lower().find(name.lower()) != -1:
                            self.atomtype[i][j] = name
        
        print('Now to change atom types')       
        while True:
            print('For each molecule, its atom types are:\n')
            count = 1
            for mol in self.atomtype:
                print('For Molecule ID: < {:d} >'.format(count))
                print(mol)
                print()
                count += 1

            print('Do you want to make any changes?, y/yes, else not')
            tget = input()
            if tget.upper() != 'Y' and tget.upper() != 'YES':
                print('Note: no changes were made...')
                print('Processing the molecules...\n')
                break
            
            print('\nWhich molecule do you want to make change?')
            print('Please input its Molecule ID')

            while True:
                tmolnm = input()
                try:
                    tmolnm = int(tmolnm)
                    if tmolnm <= 0 or len(self.atomtype) < tmolnm:
                        print('Error: no such molecule')
                        raise ValueError
                except ValueError:
                    print('Please re-input molecule ID')
                    continue
                break
                               
            print('The chosen Molecule ID is:  < {:d} >'.format(tmolnm))
            print('Its atom types are:')
            print(self.atomtype[tmolnm-1])
            print()

            print('Which atom type do you want to make a change?')
            print('Please input its Atom ID. (Counting from number 1)')
            while True:
                tatomnm = input()
                try:
                    tatomnm = int(tatomnm)
                    if tatomnm <= 0 or len(self.atomtype[tmolnm-1]) < tatomnm:
                        print('Error: no such atom')
                        raise ValueError
                except ValueError:
                    print('Please re-input atom ID')
                    continue
                break

            print('The chosen Atom is:   < {:s} >'.format(self.atomtype[tmolnm-1][tatomnm-1]))

            while True:
                print('Please input the atom name you want to use')
                tn = input()
                if len(tn) != 0: break

            print('\nThe input atom name is:    {:s}'.format(tn))
            self.atomtype[tmolnm-1][tatomnm-1] = tn
            
            print('\n\nThe final result is:\n')
            print('For Molecule ID:   < {:d} >'.format(tmolnm))
            print('Its atom types are:')
            print(self.atomtype[tmolnm-1])

            print('\n\nDo you want to continue making change? y/yes, else not')
            tp = input()
            if tp.upper() != 'Y' and tp.upper() != 'YES':
                break
            else:
                print('\n\n')



    def setorigin(self):
        """Method to rearrage atoms' sequences based on select reference"""

        # the reference list
        #
        # Format: 2D
        #
        #    [  [ nm_1,     nm_2,      nm_3,       nm_4 ], ...  ]
        #        mol_ID    Origin     X-axis      XZ-plan
        #
        #         MUST      MUST      mayHave     mayHave
        reflist = []

        print('Customize atoms\' sequences to fix atom(s) in SON calculation')
        # make a copy
        oriatomtype = self.atomtype[:]
        while True:
            if len(reflist) != 0:
                print('\n\nDo you want to continue setting? y/yes, else not')
                tp = input()
                if tp.upper() != 'Y' and tp.upper() != 'YES':
                    break
            
            print('\nFor each molecule, its atom types are:\n')
            count = 1
            for mol in self.atomtype:
                print('For Molecule ID: < {:d} >'.format(count))
                print(mol)
                print()
                count += 1
            
            print('Which molecule do you want to choose?')
            print('Please input its Molecule ID')

            while True:
                mnm = input()
                try:
                    mnm = int(mnm)
                    if mnm <= 0 or len(self.atomtype) < mnm:
                        print('Error: no such molecule')
                        raise ValueError
                except ValueError:
                    print('Please re-input molecule ID')
                    continue
                break

            print('The chosen Molecule ID is:  < {:d} >'.format(mnm))
            if len(self.atomtype[mnm-1]) == 1:
                print('Note: The chosen molecule only contains one atom')
                print('    : which will be set to Origin')
                reflist.append([mnm-1,0])
                continue

            # check reflist for repeats
            bool_repeats = False
            if mnm-1 in [i[0] for i in reflist]:
                print('Warning: for Molecule ID < {:} >, it has been already set'.format(mnm))
                print('Its setting atom types are:')
                for i,j in enumerate(self.atomtype[mnm-1]):
                    print('{:}-{:}    '.format(i+1,j),end='')
                print('\n')
                print('Do you want to override this setting? y/yes, else not')
                t = input()
                if t.upper() == 'Y' or t.upper() == 'YES':
                    bool_repeats = True
                    print('\nWarning: the old setting for Molecule ID < {:} > will be overridden\n'.format(mnm))
                    for rep in reflist:
                        if rep[0] == mnm - 1:
                            break
                    reflist.remove(rep)
                else:
                    print('Note: you have decided NOT to override old setting')
                    continue

            print('\nFor Molecule ID < {:} >'.format(mnm))
            print('Its atom types are:')
            if bool_repeats:
                for i,j in enumerate(oriatomtype[mnm-1]):
                    print('{:}-{:}    '.format(i+1,j),end='')
            else:
                for i,j in enumerate(self.atomtype[mnm-1]):
                    print('{:}-{:}    '.format(i+1,j),end='')
            print('\n')

            print('Which atoms do you want to set as reference?')
            print('Please input its Atom ID. Maximum 3 ID can be set')
            print('First-ID will be set to Origin, Second-ID will be set to positive X-axis')
            print('Third-ID will be set to XZ-plan.')
            print('Separated by one or more blank spaces. For exmpale, "1 ", " 1  2", "1  2     3"')
            while True:
                t = input()
                t = t.split()
                ls = []
                try:
                    if len(t) == 0: raise ValueError
                    if len(t) >= 4:
                        print('Error: Maximum 3 ID can be input')
                        raise ValueError
                    for i in t:
                        m = int(i)
                        if m <= 0 or len(self.atomtype[mnm-1]) < m:
                            print('Error: no such atom')
                            raise ValueError
                        ls.append(m)
                    if len(ls) != len(set(ls)):
                        print('Error: they have to be unique and exclusive to each other')
                        raise ValueError
                except ValueError:
                    print('Please re-input atom ID(s)')
                    continue
                break
            ndxlist = [mnm-1,] + [i-1 for i in ls]
            reflist.append(ndxlist)

            print('\nFor molecule ID < {:} >'.format(ndxlist[0]+1))
            if bool_repeats:
                print('Atom: < {:}-{:} > will be set to Origin'.format(ndxlist[1]+1,oriatomtype[mnm-1][ndxlist[1]]))
                if len(ndxlist) >= 3:
                    print('Atom: < {:}-{:} > will be set to positive X-axis'.format(ndxlist[2]+1,oriatomtype[mnm-1][ndxlist[2]]))
                    if len(ndxlist) >= 4:
                        print('Atom: < {:}-{:} > will be set to XZ-plan'.format(ndxlist[3]+1,oriatomtype[mnm-1][ndxlist[3]]))
            else:
                print('Atom: < {:}-{:} > will be set to Origin'.format(ndxlist[1]+1,self.atomtype[mnm-1][ndxlist[1]]))
                if len(ndxlist) >= 3:
                    print('Atom: < {:}-{:} > will be set to positive X-axis'.format(ndxlist[2]+1,self.atomtype[mnm-1][ndxlist[2]]))
                    if len(ndxlist) >= 4:
                        print('Atom: < {:}-{:} > will be set to XZ-plan'.format(ndxlist[3]+1,self.atomtype[mnm-1][ndxlist[3]]))

            self.atomtype = []
            for ndx, mol in enumerate(oriatomtype):
                bo = False
                for ndxlist in reflist:
                    if ndxlist[0] == ndx:
                        bo = True
                        break
                if bo:                        
                    ls = []
                    for t in ndxlist[1:]:
                        ls.append(mol[t])
                    for s, atom in enumerate(mol):
                        if s not in ndxlist[1:]:
                            ls.append(atom)
                    self.atomtype.append(ls)
                else:
                    self.atomtype.append(mol)
            print()
            print('Its new atom sequence is:')
            print(self.atomtype[mnm-1])
        
        # Printout final information
        print('\nAfter Atom Sequence Customization\n')
        for i, mol in enumerate(self.atomtype):
            print('For molecule ID < {:} >'.format(i+1))
            print('Its atom types are:')
            for i,j in enumerate(mol):
                print('{:}-{:}    '.format(i+1,j),end='')
            print('\n')
        print()
            
        # Now based on reflist, update self.corlist
        corlist = []
        for ndx, res in enumerate(self.corlist):
            bo = False
            for ndxlist in reflist:
                if ndxlist[0] == ndx:
                    bo = True
                    break
            if bo:
                ls = []
                for mol in res:
                    lp = []
                    for t in ndxlist[1:]:
                        lp.append(mol[t])
                    for s, atom in enumerate(mol):
                        if s not in ndxlist[1:]:
                            lp.append(atom)
                    ls.append(lp)
                corlist.append(ls)
            else:
                corlist.append(res)
        self.corlist = corlist

    

    @staticmethod
    def overlap(corlist):
        """Overlaping input coordinates to its origin by translating and rotating
           Input corlist parameter format:
               corlist: 4D
                    [     [     [    [ x, y, z ]     ]    ]    ]     
                         res   mol  atom
        """
        procorlist = []
        for res in corlist:
            if len(res[0]) == 1:
                lp = [[[0,0,0],] for mol in res]
                procorlist.append(lp)
               
            elif len(res[0]) == 2:
                lp = []
                for mol in res:
                    rx = mol[1][0] - mol[0][0]
                    ry = mol[1][1] - mol[0][1]
                    rz = mol[1][2] - mol[0][2]
                    dist12 = math.sqrt(rx*rx + ry*ry + rz*rz)             
                    lp.append( [[0,0,0],[dist12,0,0]] )
                procorlist.append(lp)
                
            elif len(res[0]) == 3:
                lp = []
                for mol in res:               
                    rx = mol[1][0] - mol[0][0]
                    ry = mol[1][1] - mol[0][1]
                    rz = mol[1][2] - mol[0][2]
                    sqr_dist12 = rx*rx + ry*ry + rz*rz
                    
                    rx = mol[2][0] - mol[0][0]
                    ry = mol[2][1] - mol[0][1]
                    rz = mol[2][2] - mol[0][2]
                    sqr_dist13 = rx*rx + ry*ry + rz*rz

                    rx = mol[2][0] - mol[1][0]
                    ry = mol[2][1] - mol[1][1]
                    rz = mol[2][2] - mol[1][2]
                    sqr_dist23 = rx*rx + ry*ry + rz*rz

                    # the base_vector coordinates, make it at the positive xz-planar
                    dist12 = math.sqrt(sqr_dist12)
                    x3 = abs( (sqr_dist13 + sqr_dist12 - sqr_dist23) / dist12 / 2 )
                    z3 = math.sqrt(sqr_dist13 - x3 * x3)

                    lp.append( [[0,0,0],[dist12,0,0],[x3,0,z3]] )
                procorlist.append(lp)

            else:
                # translate molecule's first atom to (0,0,0)
                lp1 =[]
                for mol in res:
                    ls = [[0,0,0]]
                    for atom in mol[1:]:
                        rx = atom[0] - mol[0][0]
                        ry = atom[1] - mol[0][1]
                        rz = atom[2] - mol[0][2]
                        ls.append( [rx, ry, rz] )
                    lp1.append(ls)

                # Make mol 2x >= 0. Rotate 180° around z-axis if its second atom x < 0.
                lp2 = []
                for mol in lp1:
                    if mol[1][0] < 0:
                        bool_rotate = True
                    else:
                        bool_rotate = False

                    ls = [ [0,0,0], ]
                    for atom in mol[1:]:
                        if bool_rotate:
                            rx = -atom[0]
                            ry = -atom[1]
                        else:
                            rx = atom[0]
                            ry = atom[1]
                        ls.append([rx,ry,atom[2]])
                    lp2.append(ls)

                # Make mol 2z = 0. Rotate around x-axis to make second atom to xy plane.
                # Matrix = [ [1,    0,       0  ],
                #            [0,  cos(a), sin(a)],
                #            [0, -sin(a), cos(a)]  ]
                lp3 = []
                for mol in lp2:                 
                    if mol[1][1] == 0:
                        a = math.pi/2
                    else:
                        a = math.atan2(mol[1][2], mol[1][1])

                    ry = mol[1][1]*math.cos(a) + mol[1][2]*math.sin(a)
                    ls = [[0,0,0],[mol[1][0],ry,0],]
                    
                    for atom in mol[2:]:
                        ry = atom[1]*math.cos(a) + atom[2]*math.sin(a)
                        rz = -atom[1]*math.sin(a) + atom[2]*math.cos(a)
                        ls.append([atom[0],ry,rz])
                    lp3.append(ls)

                # Make mol 2y = 0. Rotate around z-axis.
                # Matrix = [ [ cos(a),  sin(a),  0 ],
                #            [-sin(a),  cos(a),  0 ],
                #            [   0,       0,     1 ]  ]
                lp4 = []
                for mol in lp3:
                    if mol[1][0] == 0:
                        a = math.pi/2
                    else:
                        a = math.atan2(mol[1][1], mol[1][0])

                    rx = mol[1][0]*math.cos(a) + mol[1][1]*math.sin(a)
                    ls = [[0,0,0],[rx,0,0]]

                    for atom in mol[2:]:
                        rx = atom[0]*math.cos(a) + atom[1]*math.sin(a)
                        ry = -atom[0]*math.sin(a) + atom[1]*math.cos(a)
                        ls.append([rx,ry,atom[2]])
                    lp4.append(ls)


                # Make mol 3y = 0. Rotate around x-axis to make third atom to Positive xz plane.
                # Matrix = [ [1,    0,       0  ],
                #            [0,  cos(a), sin(a)],
                #            [0, -sin(a), cos(a)]  ]
                lp5 = []
                for mol in lp4:
                    if mol[2][1] == 0:
                        a = 0 if mol[2][2] >= 0 else math.pi
                    else:
                        if mol[2][2] == 0:
                            a = -math.pi/2 if mol[2][1] >= 0 else math.pi/2
                        else:
                            a = math.atan2(-mol[2][1], mol[2][2])
                    
                    rz = -mol[2][1]*math.sin(a) + mol[2][2]*math.cos(a)
                    ls = [[0,0,0],[mol[1][0],0,0],[mol[2][0],0,rz]]

                    for atom in mol[3:]:
                        ry = atom[1]*math.cos(a) + atom[2]*math.sin(a)
                        rz = -atom[1]*math.sin(a) + atom[2]*math.cos(a)
                        ls.append( [atom[0],ry,rz] )
                    lp5.append(ls)
                procorlist.append(lp5)

        return procorlist



    def file_print(self):
        if self.fout_format == 'xyz':
            self.fname = file_gen_new(self.fname,fextend='xyz')
            tot = min( [ len(i) for i in self.corlist ] )
            atomnm = sum( [ len(i[0]) for i in self.corlist ] )

            with open(self.fname,mode='wt') as f:
                for count in range(tot):
                    f.write('{:>10}\n'.format(atomnm))
                    f.write('i = {:>10}, time = {:>10}\n'.format(count,count))
                    for i in range(len(self.corlist)):
                        for j, atom in enumerate(self.corlist[i][count]):
                            f.write('{:>4} {:>14.6f} {:>14.6f} {:>14.6f}\n'.format(self.atomtype[i][j],
                                    atom[0]+self.distance*i,atom[1]+self.distance*i,atom[2]+self.distance*i))
            print('Note: generated file is < {:} >'.format(self.fname))
        elif self.fout_format == 'gro':
            pass
        elif self.fout_format == 'pdb':
            pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Molecular Spatial Overlap Normalization',allow_abbrev=False)
    parser.add_argument('-v','--version',action='version',version='mSON 0.30')
    parser.add_argument('-f','--file',help='input file path',required=True)
    parser.add_argument('--finput_format',help='Force the input file format, currently only "pdb, gro, txt" is supported')
    parser.add_argument('-d','--distance',help='Processed Molecule\'s distance, default is 5')
    parser.add_argument('-nopbc','--nopbc',help='Turn off the Periodic Boundary Condition filterion, default is False',action='store_false')
    parser.add_argument('--fixatoms',help='Custom atoms to Origin/X-axis/XZ-plan during SON calculation, default is False',action='store_true')
    parser.add_argument('--fout_format',help='Force the output file format, currently only "xyz" is supported')
    parser.add_argument('-o','--fname',help='Output file name, default is Overlap')

    args = parser.parse_args()

    fgetdict = { **vars(args) }
    bo = fgetdict.pop('nopbc',None)
    if bo:
        fgetdict['pbc'] = True
    else:
        fgetdict['pbc'] = False

    Overlap(**fgetdict).file_print()

    

