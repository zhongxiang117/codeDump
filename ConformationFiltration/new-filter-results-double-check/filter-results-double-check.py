#!/usr/bin/env python3

import os
import sys
import glob


class ReadFile:
    """
    Initialization:
        file :
            format:
                txt file
                xsf file
                gro file (under dev)
                pdb file (under dev)

    Attributes:
        system : System[Mol[Atom, ...]] : List[List[[atomtype, x,y,z], ...]]

        energy : 1D List[float]  :   None means not exist
    """
    def __init__(self,file,ext=None):
        self.file = file

        # decide file format
        if ext is None:
            ndx = file.rfind('.')
            if ndx == -1 or ndx >= len(file) - 1:
                self.ext = None
            else:
                self.ext = file[ndx+1:].lower()
            if self.ext is None: self.ext = 'txt'
        elif ext not in ['txt','xsf',]:
            print('Input file format: {:}'.format(ext))
            raise ValueError('Error: input file format is not supported')
        else:
            self.ext = ext



    def run(self):
        print('Note: reading file < {:} > ...'.format(self.file))
        self.energy = []
        if self.ext == 'txt':
            self.read_txt()
        elif self.ext == 'xsf':
            self.read_xsf()
        else:
            self.system = []



    def read_xsf(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()
        
        promol = []
        mol = []
        energylist = []
        for cnt,line in enumerate(profile):
            sub = line.strip()
            if len(sub) == 0: continue
            if sub.find('#') != -1:
                ltmp = sub.split()
                if len(ltmp) == 2:
                    energylist.append(ltmp[1])
                else:
                    energylist.append(None)
                if len(mol) == 0: continue
                promol.append(mol)
                # initialize
                mol = []
            else:
                mol.append([sub,cnt])
        # last mol
        if len(mol) != 0: promol.append(mol)
        
        if len(energylist) != len(promol):
            print('Error: energy label and coordinates are not corresponded')
            raise ValueError('Error: on file : < {:} >'.format(self.file))

        prolist = []
        for cnt,mol in enumerate(promol):
            bo = False
            if len(mol) <= 1 or mol[0][0] != 'ATOMS':
                errnum = mol[0][1] + 1
                errline = mol[0]
                bo = True
            else:
                ls = []
                for t in mol[1:]:
                    # atom info
                    atom = t[0].split()
                    if len(atom) >= 4:
                        try:
                            x = float(atom[1])
                            y = float(atom[2])
                            z = float(atom[3])
                        except ValueError:
                            bo = True
                    else:
                        bo = True
                    if bo:
                        errnum = t[1] + 1
                        errline = t[0]
                        break
                    ls.append([atom[0],x,y,z])
            if not bo:
                if energylist[cnt] is None:
                    self.energy.append(None)
                else:
                    try:
                        self.energy.append(float(energylist[cnt]))
                    except ValueError:
                        print('\nWarning: for molecule')
                        for m in ls: print(m)
                        print('Warning: energy is incorrect, ignoring\n')
                        self.energy.append(None)
            if bo:
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...')
            else:
                prolist.append(ls)
        
        self.system = self.check_atomtype(prolist)



    def read_txt(self):
        """read coordinates from txt file"""
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        # List[List[[atomtype, x, y, z], ...]]
        promol = []
        mol = []
        energylist = []
        for cnt,line in enumerate(profile):
            sub = line.strip()
            if len(sub) == 0:
                if len(mol) == 0: continue
                promol.append(mol)
                # initialize
                mol = []
            else:
                mol.append([sub,cnt])
        # last mol
        if len(mol) != 0: promol.append(mol)

        prolist = []
        for mol in promol:
            # check whether energy exist or not
            firstlabel = mol[0][0]
            ene = None
            if firstlabel[0] == '#':
                mol = mol[1:]
                ltmp = firstlabel.split()
                if len(ltmp) == 2:
                    try:
                        ene = float(ltmp[1])
                    except ValueError:
                        pass
            self.energy.append(ene)

            bo = False
            ls = []
            for t in mol:
                # atom info
                atom = t[0].split()
                if len(atom) == 4:
                    try:
                        x = float(atom[1])
                        y = float(atom[2])
                        z = float(atom[3])
                    except ValueError:
                        bo = True
                else:
                    bo = True
                if bo:
                    errnum = t[1] + 1
                    errline = t[0]
                    break
                ls.append([atom[0],x,y,z])
            if bo:
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...')
            else:
                prolist.append(ls)

        self.system = self.check_atomtype(prolist)



    def check_atomtype(self,syslist):
        """check atomtype for syslist & return results basis on first entry

        Inputs:
            syslist : System[Mol[Atoms]] : List[List[[atomtype, x, y, z], ...]]

        Return:
            same format with inputs
        """
        if len(syslist) <= 1: return syslist

        ndxlist = [i[0] for i in syslist[0]]
        prolist = [syslist[0],]
        for mol in syslist[1:]:
            bo = True
            for cnt,atom in enumerate(mol):
                if ndxlist[cnt] != atom[0]:
                    print('Wrong: not cooresponded: {:}'.format(atom))
                    print('Note: whole sets are omitted...\n')
                    bo = False
                    break
                if not bo: break
            if bo:
                prolist.append(mol)

        return prolist



USAGE = """
filter-results-double-check.py   SystemFile   -n IndexFiles
"""

if len(sys.argv) <= 1:
    print(USAGE)
    exit()

line = ' '.join(sys.argv[1:])
ltmp = line.split('-n')
if len(ltmp) != 2:
    print('Wrong: input: < {:} >'.format(line))
    exit()

sysfile = ltmp[0].split()
bo = True
if len(sysfile) != 1:
    print('Warning: wrong no input system file')
    bo = False
else:
    if os.path.isfile(sysfile[0]):
        sysfile = sysfile[0]
        ndxfile = []
        for file in ltmp[1].split():
            if os.path.isfile(file):
                ndxfile.append(file)
            else:
                print('Warning: not file < {:} >, ignoring'.format(file))
                break
        if len(ndxfile) == 0:
            print('Warning: no index files')
        else:
            bo = False
    else:
        print('Warning: not an input file < {:} >'.format(sysfile[0]))
if bo: exit()


FD = ReadFile(sysfile)
FD.run()

root_system = FD.system
root_energy = FD.energy

comp_system = []
comp_energy = []
for f in ndxfile:
    FD = ReadFile(f)
    FD.run()
    for s in FD.system: comp_system.append(s)
    for e in FD.energy: comp_energy.append(e)



for cnt,rsys in enumerate(root_system):
    rene = root_energy[cnt]
    bo = True
    for ndx,psys in enumerate(comp_system):
        bo = True
        if len(rsys) == len(psys):
            pene = comp_energy[ndx]
            if rene is None:
                if pene is not None:
                    bo = False
            else:
                if pene is None:
                    bo = False
                elif abs(rene-pene) > 0.000001:
                    bo = False
            if bo:
                for n,ra in enumerate(rsys):
                    pa = psys[n]
                    if len(ra) != len(pa) or len(ra) != 4:
                        print('Fatal: how could it be? -- atom-nms')
                        exit()
                    if pa[0] != ra[0]:
                        print('Fatal: how could it be? -- atom-label')
                        exit()
                    dx = abs(pa[1]-ra[1])
                    if dx > 0.000001:
                        bo = False
                    else:
                        dy = abs(pa[2]-ra[2])
                        if dy > 0.000001:
                            bo = False
                        else:
                            dz = abs(pa[3]-ra[3])
                            if dz > 0.000001: bo = False
        else:
            print('Fatal: how could it be? -- molnms')
            exit()
        if bo: break
    if not bo:
        print('Fatal: not found')
        print('Fatal: energy -- {:}'.format(rene))
        print('Fatal: mol --')
        for i in rsys: print(i)


print('Note: done comparing for file < {:} >'.format(sysfile))


