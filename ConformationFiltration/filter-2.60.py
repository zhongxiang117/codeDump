#!/usr/bin/env python3

import os
import sys
import math
import argparse
import matplotlib.pyplot as plt


FEATURES = [
    'version 0.10 : start',
    'version 0.20 : add bond connection',
    'version 0.30 : add angles, work on square cos',
    'version 0.40 : change angle to sin',
    'version 0.50 : split input as fragments, same on bonds & angles',
    'version 0.60 : change angle to itself, unit on degree',
    'version 0.70 : add more type for bond perception'
    'version 0.80 : refine filter function',
    'version 0.90 : codes execution efficiency',
    'version 1.00 : printout more info, RELEASE',
    'version 1.10 : add new class CrossFiltration',
    'version 1.20 : make class Filtration more robust & general',
    'version 1.30 : add filtration results save to image',
    'version 1.40 : add bulk process',
    'version 1.50 : bulk process is perfectly done, RELEASE',
    'version 1.60 : add argparse, RELEASE',
    'version 1.70 : printout more info, RELEASE',
    'version 1.80 : add argparse features, RELEASE',
    'version 1.90 : save probability bins',
    'version 2.00 : add Boolean for images output',
    'version 2.10 : add supprot for Backup xsf file',
    'version 2.20 : fix ReadFile to read the last molecule',
    'version 2.21 : fix error when check user input connections',
    'version 2.22 : print out final length & ratio for BulkProcess',
    'version 2.30 : preserve energy for xsf file',
    'version 2.31 : fix BulkProcess fname',
    'version 2.40 : add energy for ReadFile read_txt',
    'version 2.50 : add SaveFile xsf',
    'version 2.60 : add file format explanation',
]


VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


FILEFORMAT = """
Input File Format for BOSS Output Filtration

Currently, only two types of file formats are supported.

For txt file:

    molecules are separated by space(s), char '#' can be omitted

            [   #   energy
            |   atom-1   x    y    z
    mol     |   atom-2   x    y    z
            |   atom-3   x    y    z
            |   ...
            [   <space>


For xsf file:

    molecules are separated by '#', keyword 'ATOMS' is important,
    case sensitive, otherwise, the whole set will be ignored,
    follwing their warning info, to be used for double checking

            [   #   energy
            |
            |   ATOMS
    mol     |   atom-1   x    y    z    else
            |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...
"""


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




class BondPerception:
    """Bond connecting perception

    Paper:
        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Note:
        1) only the first Detect Connection Rule is used
        2) limited for atom H, C, N, O, F, S, Cl

    Inputs:
        mol : single : Mol[[type,x,y,z], ...] : List[[type,x,y,z], ...]

        refer : which bond type     : default bonded, one of [1,'bonded',True]
            int:0 | str:all         | None  :   choose all atom bond connection
            int:1 | str:bonded      | True  :   choose only bonded
            int:2 | str:nonbonded   | False :   choose only non-bonded
            int:3 | str: fnb,fragments,piece :   nonbonded based on fragments

    Attributes:
        bcon        : System[[atom-i, atom-j], ...] : List[[int,int], ...]

        fragments   : valid when refer=3    : 2D : List[ List[int], ...]
    """
    RADIUS = {
        'H' : 0.23,
        'C' : 0.68,
        'N' : 0.68,
        'O' : 0.68,
        'F' : 0.64,
        'S' : 1.02,
        'Cl': 0.99,
    }
    ATOMTYPE = {
        '1' : 'H',
        '6' : 'C',
        '7' : 'N',
        '8' : 'O',
        '9' : 'F',
        '16': 'S',
        '17': 'Cl',
    }
    def __init__(self,mol,refer=True):
        self.mol = mol
        if refer in [0,'all',None]:
            self.refer = 0
        elif refer in [1,'bonded',True]:
            self.refer = 1
        elif refer in [2,'nonbonded',False]:
            self.refer = 2
        elif refer in [3,'fnb','fragments','piece']:
            self.refer = 3
        else:
            raise ValueError('not defined: refer')



    def run(self):
        """calc bond index based on self.refer"""
        print('Note: calculating bond perception ...')
        self.bcon = []
        self.fragments = []
        if len(self.mol) <= 1: return

        fall = []
        for i,tmp in enumerate(self.mol[:-1]):
            j = i + 1
            while j < len(self.mol):
                fall.append([i,j])
                j += 1

        if self.refer == 0:
            self.bcon = fall
            return

        dislist = self.calc_input_atom_square_distance()
        refdict = self.calc_refer_atom_square_distance()

        fbond = []
        for ndx in dislist:
            ref = refdict[ndx[0]]
            if ndx[1] >= ref[0] and ndx[1] <= ref[1]:
                t = ndx[0].split('-')
                fbond.append([int(t[0]),int(t[1])])

        if self.refer == 1:
            self.bcon = fbond
            return

        fnon = [i for i in fall if i not in fbond]
        if self.refer == 2:
            self.bcon = fnon
            return

        if len(fbond) == 0:
            self.bcon = fall
            return

        fdict = {}
        for k in range(len(self.mol)):
            ls = [i[1] for i in fbond if i[0] == k]
            if len(ls) != 0:
                fdict[k] = ls

        fref = []
        for k in range(len(self.mol)):
            if k in fdict:
                ls = fdict[k]
                ls.append(k)
                fdict.pop(k)
                while True:
                    ref = [t for t in ls if t in fdict]
                    if len(ref) == 0:
                        break
                    for i in ref:
                        # list calculation
                        ls += fdict[i]
                        fdict.pop(i)
                fref.append(sorted(ls))

        fndx = []
        for ndx in fref:
            for i,k in enumerate(ndx[:-1]):
                j = i + 1
                while j < len(ndx):
                    fndx.append([k,ndx[j]])
                    j += 1

        fnon3 = [i for i in fall if i not in fndx]
        if self.refer == 3:
            self.bcon = fnon3
            self.fragments = fref
            return



    def calc_input_atom_square_distance(self):
        """
        Rule:
            calculation is performed based on sequence (starting at 0)

            For example, if input list indexed as [0, 1, 2, 3],
            Then square distance [0-1, 0-2, 0-3, 1-2, 1-3, 2-3] is calculated

        Return:
            dislist : [[index-index, SD], ...] : List[[str,float], ...]
        """
        if len(self.mol) <= 1: return []

        dislist = []
        for i, ref in enumerate(self.mol[:-1]):
            j = i + 1
            while j < len(self.mol):
                atom = self.mol[j]
                dx = ref[1] - atom[1]
                dy = ref[2] - atom[2]
                dz = ref[3] - atom[3]
                ds = dx*dx + dy*dy + dz*dz
                ndx = '{:}-{:}'.format(i,j)
                dislist.append([ndx,ds])
                j += 1

        return dislist



    def calc_refer_atom_square_distance(self):
        """
        Rule:
            calculation is performed based on sequence (starting at 0)

            For example, if input list indexed as [0, 1, 2, 3],
            Then square distance [0-1, 0-2, 0-3, 1-2, 1-3, 2-3] is calculated

        Return:
            refdict : { index-index: [SDlow, SDhigh], ...} : Dict{str:List, ...}

        Equation:
            0.8 < dij < ri + rj + 0.4
        """
        if len(self.mol) <= 1: return {}

        radiuslist = []
        bo = False
        for atom in self.mol:
            if isinstance(atom[0],int):
                if str(atom[0]) in self.ATOMTYPE:
                    at = self.RADIUS[self.ATOMTYPE[str(atom[0])]]
                else:
                    bo = True
            elif isinstance(atom[0],str):
                if atom[0] in self.RADIUS:
                    at = self.RADIUS[atom[0]]
                elif str(atom[0]) in self.ATOMTYPE:
                    at = self.RADIUS[self.ATOMTYPE[str(atom[0])]]
                else:
                    bo = True
            else:
                bo = True
            if bo:
                print('Error: atom type not defined: {:}'.format(atom[0]))
                raise ValueError('not defined: atom type')
            radiuslist.append(at)

        refdict = {}
        for i, ref in enumerate(radiuslist[:-1]):
            j = i + 1
            while j < len(radiuslist):
                high = ref + radiuslist[j] + 0.4
                ndx = '{:}-{:}'.format(i,j)
                refdict[ndx] = [0.64, high*high]
                j += 1

        return refdict




class Filtration:
    """Filter molecules based on Bond-Connection

    Inputs:
        system : System[Mol[Atoms]] : List[List[[atomtype, x, y, z], ...]]
        bcon : System[[atom-i, atom-j], ...] : List[[int,int], ...]
        acon : System[[atom-i, atom-j, atom-k], ...] : List[[int,int,int], ...]

        obonds  :  Boolean  :  whether calculate bonds  :  default False
        oangles :  Boolean  :  whether calculate angles :  default False
        =>  opar  :  Boolean  :  whether calculate prob_par  :  default False
        =>  oall  :  Boolean  :  whether calculate prob_all  :  default True


    Attributes:
        system  :  good molecules after filtration
        sysbad  :  filtered out molecules

        prob_initial : initial probability  :   dict

            keys:
                bonds   :   List
                    prob_par    : 3D : List[ [List[int],float] ]
                    prob_all    : 2D : List[ List[int],  float ]

                angles  :   List
                    prob_par    : 3D : List[ [List[int],float] ]
                    prob_all    : 2D : List[ List[int],  float ]

        prob_final :  same format as prob_initial

        tolerance  :  [bond, angle]  :  [Angstrom, Degree]  : default [0.1,0.1]

        bondlist   :  2D  :  List[ List[float] ]  :  good, correspond to bcon
        anglelist  :  2D  :  List[ List[float] ]  :  good, correspond to acon

        reflist    :  1D  List[int]     :   index of for bad molecules


    Caution:
            Because the BOSS solute move is only for the single time move,
            so Sum(difference for untouched atoms) =~ 0.0,
            thus, when doing filtration, please take care of <dt>
    """
    def __init__(self,system,*args,**kwargs):
        self.system = system

        self.bcon = kwargs['bcon'] if 'bcon' in kwargs else None
        self.acon = kwargs['acon'] if 'acon' in kwargs else None
        self.bondlist = []
        self.anglelist = []

        if 'tolerance' in kwargs and kwargs['tolerance'] is not None:
            # Angstrom, Degree
            self.tolerance = kwargs['tolerance']
        else:
            self.tolerance = [0.1, 0.1]

        obonds = kwargs['obonds'] if 'obonds' in kwargs else None
        oangles = kwargs['oangles'] if 'oangles' in kwargs else None
        opar = kwargs['opar'] if 'opar' in kwargs else None
        oall = kwargs['oall'] if 'oall' in kwargs else None

        if obonds is None:
            self.obonds = False if self.bcon is None else True
        else:
            self.obonds = True if obonds is True else False

        if oangles is None:
            self.oangles = False if self.acon is None else True
        else:
            self.oangles = True if oangles is True else False

        self.opar = True if opar is True else False
        self.oall = True if oall is True else False

        self.prob_initial = {'angles':[], 'bonds':[]}
        self.prob_final = {'angles':[], 'bonds':[]}
        self.sysbad = []



    def run(self):
        """attemption on filtering

        Attributes:
            self.system : update
            self.sysbad : filter out molecules
        """
        # to improve efficiency, bondlist only needs to be calculated once
        bondlist = self.calc_square_bond_connection_distance()
        anglelist = self.calc_angle_degree_values()

        # increments
        binc = self.tolerance[0] * self.tolerance[0]
        ainc = self.tolerance[1]

        # initial
        if self.obonds or self.oangles:
            print('Note: calculating initial probability ...')
            self.probability(self.prob_initial,bondlist,anglelist,binc,ainc)

        reflist = self.calc_mols_filter_list(bondlist,anglelist,[binc,ainc])
        self.reflist = reflist

        # update self.system
        # => equals to update bondlist & anglelist
        #
        # Similarly like:
        #   if we have a list: [a, b, c, d, e]
        #   we want to remove values at index [1,2,4], e.g. [b,c,e]
        #   how can we do that?
        #
        # to increase efficiency, switch system, sysbad accordingly.
        print('Note: updating system ...')
        self.sysbad = [ndx for i,ndx in enumerate(self.system) if i in reflist]
        self.system = [ndx for i,ndx in enumerate(self.system) if i not in reflist]

        print('Note: updating bondlist & anglelist ...')
        self.bondlist = [ndx for i,ndx in enumerate(bondlist) if i not in reflist]
        self.anglelist = [ndx for i,ndx in enumerate(anglelist) if i not in reflist]

        # final
        if self.obonds or self.oangles:
            print('Note: calculating final probability')
            self.probability(self.prob_final,self.bondlist,self.anglelist,binc,ainc)


    def probability(self,fdict,bondlist,anglelist,binc,ainc):
        """probability exclusively for self.prob_initial & self.prob_final"""
        if self.obonds:
            print('  ==> bonds ...')
            fdict['bonds'] = self.calc_probability_bonds(bondlist,binc)

        if self.oangles:
            print('  ==> angles ...')
            fdict['angles'] = self.calc_probability_angles(anglelist,ainc)



    def calc_probability_angles(self,anglelist,dt):
        """
        Return:
            prob_par  : 3D : List[ [List[int],float] ]
            prob_all  : 2D : List[ List[int],  float ]

        Similar rules like bonds, please refer it for more detail
        """
        return self.calc_probability_bonds(anglelist,dt)



    def calc_probability_bonds(self,bondlist,dt):
        """
        Inputs:
            dt : float

        Return:

            prob_par  : 3D : List[ [List[int],float] ]

                => means for any corresponding bond-length in Bond-Connection,
                   the probability of any two molecules whose difference fall
                   within defined range

                explanation:
                    assume three molecules, the number of atoms are 4;
                    A = [a1, a2, a3, a4]
                    B = [b1, b2, b3, b4]
                    C = [c1, c2, c3, b4]

                    define its Bond-Connection index is: [[0,1], [1,2], [1,3]]
                    means bonded atom pairs in A are [a1-a2, a2-a3, a2-a4]

                    now, calculate their bondlist;

                    DA = [Da01, Da12, Da13]
                    DB = [Db01, Db12, Db13]
                    DC = [Dc01, Dc12, Dc13]

                    thus, the difference of bond-length for any two molecules
                    corresponding to Bond-Connection will be DA-DB, DA-DC, DB-DC

                    therefore a new array will be got:
                    New = [Da01-Db01, Da01-Dc01, Db01-Dc01]

                    finally we can calculate molecule-bond-length probability


            prob_all  : 2D : List[ List[int],  float ]

                => means probability on overall bonds difference

                explanation:
                    continuing top, we have already calculated DA, DB, DC

                    calculate their average values;

                    TA = Average(DA)
                    TB = Average(DB)
                    TC = Average(DC)

                    following the same rule, calculate differences,
                    New = [TA-TB, TA-TC, TB-TC]

                    finally we can calculate overall molecule probability


        Note:
            1) calculation is based on square values
            2) not normalized
            3) List[float] is corresponded with self.bcon
        """
        def calc_list(ls,dt):
            rmin = min(ls)
            steps = (max(ls)- rmin) / dt
            steps = int(steps+0.5)

            # all on same values
            if steps == 0: return [len(ls)],rmin

            prolist = []
            for v in range(steps):
                low = rmin + v*dt
                high = rmin + (v+1)*dt
                cnt = 0
                for i in ls:
                    if i >= low and i < high:
                        cnt += 1
                prolist.append(cnt)

            return prolist,rmin


        if len(bondlist) <= 1: return [],[]

        # calc individual probability for each bonds
        prob_par = []
        if self.opar:
            print('    --> computation on par entry ...')
            for cnt in range(len(bondlist[0])):
                ls = [t[cnt] for t in bondlist]
                p = calc_list(ls,dt)
                prob_par.append(p)

        prob_all = []
        if self.oall:
            print('    --> computation on molecules ...')
            # Caution! only one movement
            #sub = [sum(i)/len(i) for i in bondlist]
            sub = [sum(i) for i in bondlist]
            prob_all = calc_list(sub,dt*len(bondlist[0]))

        return prob_par, prob_all



    def calc_mols_filter_list(self,bondlist,anglelist,tolerance):
        """
        Inputs:
            tolerance  :  [bond^2, angle]  :  List

        Rule:
            Since the variance is only based on the single movement,
            either can be a bond or angle, so the changes on any two adjacent
            molecules are larger than tolerance will be removed

        Caution:
            Because the BOSS solute move is only for the single time move,
            thus, dt should not multiple the total numbers,
            because Sum(Amol untouched atoms) =~ 0.0

        Return:
            reflist  :  List[int]  :  index of molecules waiting to be removed
        """
        # note: bondlist or anglelist may not exist
        # always in sequence: bl, al
        if len(bondlist) != 0:
            bl = [sum(i) for i in bondlist]
            if len(anglelist) != 0:
                al = [sum(i) for i in anglelist]
                if len(bl) != len(al):
                    print('How can it be?')
                    raise ValueError('not correspond, Weird')
            else:
                al = []
        else:
            if len(anglelist) != 0:
                bl = [sum(i) for i in anglelist]
            else:
                bl = []
            al = []

        if len(bl) <= 3: return []

        print('Note: calculating filtering list ...')
        # sort by index
        ndxlist = sorted(range(len(bl)),key=lambda k: bl[k])

        # deal with caution 2
        odd1 = ndxlist[-1]
        odd2 = ndxlist[-2]

        # This part does is to recursively remove index, whose difference
        # separately with two adjacent values is smaller than tolerance
        #
        # for example, we have inputs like;
        #
        #   bl      = [0.0, 0.1, 0.15, 0.18, 0.19, 0.20, 0.3, 0.4, 0.43, 0.5]
        #   al      = [0.0, 0.1, 0.15, 0.18, 0.19, 0.20, 0.3, 0.4, 0.43, 0.5]
        # ndxlist   =   0    1    2     3     4     5     6    7    8     9
        #
        # tolerance = [0.1, 0.1]
        #
        # then we will know if we remove index in [2,3,4,8], the new array
        # blnew = [0.0, 0.1, 0.20, 0.3, 0.4, 0.5] will meet the requirement
        #
        # Caution:
        #   1) round off error is ignored
        #   2) if the last two are bad, they will be both removed
        reflist = []
        while True:
            bdel = [bl[j] - bl[ndxlist[i]] for i,j in enumerate(ndxlist[1:])]
            ls = []
            cnt = 0
            while cnt < len(bdel):
                bo = False
                if bdel[cnt] < tolerance[0]:
                    ndx = cnt + 1

                    # care when al not exist
                    if len(al) != 0:
                        da = al[ndx] - al[cnt]
                        if da < tolerance[1]:
                            bo = True
                        else:
                            cnt += 1
                    else:
                        bo = True
                if bo:
                    ls.append(ndx)
                    cnt += 2
                else:
                    cnt += 1

            if len(ls) == 0:
                break

            for i in ls: reflist.append(ndxlist[i])

            ndxlist = [ndx for i,ndx in enumerate(ndxlist) if i not in ls]
            if len(ndxlist) <= 1:
                break

        if odd1 in reflist and odd2 in reflist: reflist.remove(odd1)

        return reflist



    def calc_square_bond_connection_distance(self):
        """calc atom pairs distance square based on Bond-Connection

        Return:
            bondlist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]

            [l1, l2, ...] : distance square value defined in self.bcon
        """
        print('Note: calculating bondlist ...')
        if len(self.system) == 0 or len(self.system[0]) <= 1: return []

        # note: in square value
        bondlist = []
        if self.bcon is None:
            for mol in self.system:
                at1 = mol[0]
                ls = []
                for at2 in mol[1:]:
                    dx = at1[1] - at2[1]
                    dy = at1[2] - at2[2]
                    dz = at1[3] - at2[3]
                    tmp = dx*dx + dy*dy + dz*dz
                    ls.append(tmp)
                bondlist.append(ls)

            self.bcon = [[0,i+1] for i in range(len(self.system[0])-1)]
        elif len(self.bcon) != 0:
            for mol in self.system:
                ls = []
                for ndx in self.bcon:
                    at1 = mol[ndx[0]]
                    at2 = mol[ndx[1]]
                    dx = at1[1] - at2[1]
                    dy = at1[2] - at2[2]
                    dz = at1[3] - at2[3]
                    tmp = dx*dx + dy*dy + dz*dz
                    ls.append(tmp)
                bondlist.append(ls)

        return bondlist



    def calc_angle_degree_values(self):
        """Angle list

        Note:
            Based on self.acon

            Or;

            It is based on the first input first two atoms as the reference.
            For example, if we have a system contains many atoms bigger than 2,
            sequentially, denote them as; A (1st), B (2nd), E (else)...

            Then, angles, E-A-B will be calculated. Thus starting from third
            atom, we have an new array, [A3, A4, A5, ...]

        Rule:
            assume cooridnates,

            A (ax, ay, az)
            B (bx, by, bz)
            E (ex, ey, ez)

            Vector,

            AB = (ax-bx, ay-by, az-bz)
            AE = (ax-ex, ay-ey, az-ez)

            cos<EAB> = Sum(abi*aei) / dis(AB) * dis(AE)

            <EAB> = math.acos(sigma) * 180.0 / math.pi

        Return:
            anglelist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]

            empty list for number of atoms are less than 3
        """
        print('Note: calculating anglelist ...')
        if len(self.system) == 0 or len(self.system[0]) <= 2: return []

        # note: in square value
        anglelist = []
        cvt = 180.0 / math.pi
        if self.acon is None:
            for mol in self.system:
                A = mol[0]
                B = mol[1]
                AB = [A[1]-B[1],A[2]-B[2],A[3]-B[3]]
                ls = []
                for E in mol[2:]:
                    AE = [A[1]-E[1],A[2]-E[2],A[3]-E[3]]
                    tot = AB[0]*AE[0] + AB[1]*AE[1] + AB[2]*AE[2]
                    sub = sum([i*i for i in AB]) * sum([i*i for i in AE])
                    rst = math.acos(tot/pow(sub,0.5)) * cvt
                    ls.append(rst)
                anglelist.append(ls)

            self.acon = [[0,1,i+2] for i in range(len(self.system[0])-2)]
        elif len(self.acon) != 0:
            for mol in self.system:
                ls = []
                for ndx in self.acon:
                    A = mol[ndx[0]]
                    B = mol[ndx[1]]
                    E = mol[ndx[2]]
                    AB = [A[1]-B[1],A[2]-B[2],A[3]-B[3]]
                    AE = [A[1]-E[1],A[2]-E[2],A[3]-E[3]]
                    tot = AB[0]*AE[0] + AB[1]*AE[1] + AB[2]*AE[2]
                    sub = sum([i*i for i in AB]) * sum([i*i for i in AE])
                    rst = math.acos(tot/pow(sub,0.5)) * cvt
                    ls.append(rst)
                anglelist.append(ls)

        return anglelist




def func_calc_connection(fragments):
    """Bond & Angle connection based on fragments

    Input:
        fragments  :  2D  :  List[ List[int] ]

    Note:
        Based on first 2 atoms, which is chosen by sequence,
        find bond & angle connection with all atoms in other fragments

    Return:
        bcon  :  2D  :  List[ List[int,int], ...]
        acon  :  2D  :  List[ List[int,int,int], ...]
    """
    print('Note: calculating bond & angle connections ...')
    if len(fragments) <= 1: return [],[]

    bcon = []
    b1 = fragments[0][0]
    for ndx in fragments[1:]:
        for v in ndx:
            bcon.append([b1,v])

    acon = []
    for ndx,mol in enumerate(fragments):
        if len(mol) >= 2:
            break
    if len(mol) >= 2:
        at1 = mol[0]
        at2 = mol[1]
        for cnt,mol in enumerate(fragments):
            if cnt != ndx:
                for v in mol:
                    acon.append([at1,at2,v])
    else:
        if len(fragments) >= 3:
            at1 = fragments[0][0]
            at2 = fragments[1][0]
            for v in fragments[2:]:
                acon.append([at1,at2,v[0]])

    return bcon,acon




class CrossFiltration(Filtration):
    """filter out system based on sysndx

    Inputs:
        system : System[Mol[Atoms]] : List[List[[atomtype, x, y, z], ...]]
        sysndx : System[Mol[Atoms]] : List[List[[atomtype, x, y, z], ...]]
        tolerance : 1D : List[float, float] : [bond, angle]

    Explanation:
        there are two systems, we want to use sysndx as the index filter out
        all similar molecule in system, based on given tolerance.

    Note:
        if sysndx is optional, self filtration will always be on
    """
    def __init__(self,system,sysndx=None,*args,**kwargs):
        super().__init__(system,*args,**kwargs)
        self.sysndx = sysndx

    def run(self):
        """overwrite parent method"""
        # turn off all probability calculation
        self.obonds = False
        self.oangles = False
        self.opar = False
        self.oall = False

        # self filtration
        super().run()

        if self.sysndx is None:
            # filtration is done
            # but we still initialize sysndx to empty list
            self.sysndx = []
        else:
            # cross filtration

            # shallow copy
            sysgood = self.system

            # 1st, for sysndx, calculate bondlist & anglelist
            ndxlth = len(self.sysndx)
            self.system = self.sysndx
            bltot = self.calc_square_bond_connection_distance()
            altot = self.calc_angle_degree_values()

            # 2nd, calculation reflist
            binc = self.tolerance[0] * self.tolerance[0]
            ainc = self.tolerance[1]

            # append data onto sysndx results
            for i in self.bondlist: bltot.append(i)
            for i in self.anglelist: altot.append(i)

            # 3rd, calculate reflist
            reflist = self.calc_mols_filter_list(bltot,altot,[binc,ainc])
            reflist = [i-ndxlth for i in reflist if i > ndxlth]

            # 3rd, update
            self.system = []
            for i,mol in enumerate(sysgood):
                if i in reflist:
                    self.sysbad.append(mol)
                else:
                    self.system.append(mol)

            self.bondlist = [mol for i,mol in enumerate(self.bondlist) if i not in reflist]
            self.anglelist = [mol for i,mol in enumerate(self.anglelist) if i not in reflist]




def file_gen_new(fname,fextend='txt',foriginal=True,bool_dot=True):
    """Generate new file name without overwritings

    Args:
        fname   (str)   :   input file fname
        fextend (str)   :   file extension
        foriginal (bool):   whether keep original
        bool_dot (bool) :   force check dot convention or not

    Returns:
        str     :   new file name
    """
    filename = fname
    pos = filename.rfind('.')
    if bool_dot and pos != -1:
        fname = filename[:pos]
        fextend = filename[pos:]
    else:
        fextend = '.' + fextend

    if foriginal is True:
        if not os.path.isfile(fname+fextend):
            return fname+fextend

    i = 1
    filename = fname
    while True:
        fname = filename + '-' + str(i) + fextend
        if not os.path.isfile(fname): break
        i += 1
    return fname




def plot_filtration_save_image(ini,fin=None,dt=None,fname=None,key=None):
    """
    Inputs:
        ini     :   2D  :   List [ List[int, ...],  float]
        fin     :   2D  :   List [ List[int, ...],  float]
        dt      :   float   :   increments, optional
        fname   :   str :   warning, overwritten may happen
        key     :   str :   {bonds, angles}, specify plot type

    Return:
        True if file is successfully generated, otherwise, False
    """
    if fname is None: fname = 'filtration-image.jpg'

    boi = True
    if ini is None or len(ini) == 0: boi = False

    bof = True
    if fin is None or len(fin) == 0: bof = False

    if (not boi) and (not bof): return False

    # check data size
    if len(ini[0]) <= 5: boi = False
    if len(fin[0]) <= 5: bof = False
    if (not boi) and (not bof): return False

    if boi:
        yini = ini[0]
    if bof:
        yfin = fin[0]

    if dt is None:
        if boi:
            xini = list(range(len(yini)))
        if bof:
            xfin = list(range(len(yfin)))
    else:
        if boi:
            xini = [ini[1]+dt*t for t in range(len(yini))]
        if bof:
            xfin = [fin[1]+dt*t for t in range(len(yfin))]

    if key is None:
        title = 'Filtration'
    else:
        if dt is None:
            if key.lower() == 'bonds':
                title = 'Filtration on Bond (Angstrom)'
            elif key.lower() == 'angles':
                title = 'Filtration on Angle (Degree)'
            else:
                title = 'Filtration'
        else:
            if key.lower() == 'bonds':
                title = 'Filtration on Bond ({:} Angstrom)'.format(dt)
            elif key.lower() == 'angles':
                title = 'Filtration on Angle ({:} Degree)'.format(dt)
            else:
                title = 'Filtration'

    if boi and bof:
        # both exist
        plt.plot(xini,yini,'r-',xfin,yfin,'b-')
    elif boi:
        plt.plot(xini,yini,'r-')
    elif bof:
        plt.plot(xfin,yfin,'r-')
    else:
        print('Fatal: this should be never executed')
    plt.title(title)

    # now save figure
    print('Note: filtration plot is saved to file < {:} >'.format(fname))
    plt.savefig(fname)
    plt.close()
    return True




class SaveFile:
    """opposite operation to ReadFile

    Inputs:
        system : 2D  :  Mol[Atom, ...]  :  List[[atomtype, x,y,z], ...]

        ftype  : Output file type
            format
                txt file
                gro file    (under dev)
                pdb file    (under dev)

        force_real_atom_type     :   Boolean    :   default False

        fname : file to be saved, warning, overwritten may happen
    """
    def __init__(self,system,*args,**kwargs):
        if len(system) == 0:
            print('Warning: no inputs')
            return
        self.system = system

        if 'ftype' in kwargs and kwargs['ftype'] is not None:
            self.ftype = kwargs['ftype'].lower()
        else:
            self.ftype = 'txt'
        if self.ftype not in ['txt','xsf']:
            print('Warning: currently not support < {:} >'.format(self.ftype))
            raise ValueError('not support')

        if 'force_real_atom_type' in kwargs:
            ft = kwargs['force_real_atom_type']
        else:
            ft = False
        self.force_real_atom_type = True if ft is True else False

        if 'fname' in kwargs and kwargs['fname'] is not None:
            if len(kwargs['fname'].split()) != 0:
                self.fname = kwargs['fname']
            else:
                self.fname = 'system.txt'
        else:
            self.fname = 'system.txt'

        if 'energy' in kwargs and kwargs['energy'] is not None:
            self.energy = kwargs['energy']
        else:
            self.energy = []



    def run(self):
        """save to file"""
        atlist = self.get_atomtype()
        if self.ftype == 'txt':
            self.save_txt(self.system,atlist,self.fname)
        elif self.ftype == 'xsf':
            self.save_xsf(self.system,atlist,self.fname)



    def save_xsf(self,system,atlist,fname):
        print('Note: saving to file < {:} >'.format(fname))
        with open(fname,'wt') as f:
            for ndx,mol in enumerate(system):
                if len(self.energy) != 0 and self.energy[ndx] is not None:
                    line = '#  {:}\n'.format(self.energy[ndx])
                else:
                    line = '#\n'
                f.write(line)
                f.write('\n')
                f.write('ATOMS\n')
                ndx = 0
                for i in mol:
                    line = '{:3} {:>10} {:>10} {:>10}   1.0  1.0  1.0\n'.format(atlist[ndx],i[1],i[2],i[3])
                    f.write(line)
                    ndx += 1
                f.write('\n\n')



    def save_txt(self,system,atlist,fname):
        print('Note: saving to file < {:} >'.format(fname))
        with open(fname,'wt') as f:
            for ndx,mol in enumerate(system):
                if len(self.energy) != 0 and self.energy[ndx] is not None:
                    f.write('#  {:}\n'.format(self.energy[ndx]))

                for cnt,at in enumerate(mol):
                    line = '{:<2} {:>15} {:>15} {:>15}\n'.format(atlist[cnt],*at[1:])
                    f.write(line)
                f.write('\n\n')



    def get_atomtype(self):
        pt_nm = {
            1:'H', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 15:'P', 16:'S', 17:'Cl',
            35:'Br', 53:'I',
        }

        pt_ch = {}
        for k,v in pt_nm.items(): pt_ch[str(k)] = v

        atlist = [j[0] for i in self.system for j in i]
        if self.force_real_atom_type:
            print('Note: getting real atom type ...')
            reflist = []
            for at in atlist:
                if at in pt_nm:
                    at = pt_nm[at]
                elif at in pt_ch:
                    at = pt_ch[at]
                reflist.append(at)
            atlist = reflist

        return atlist




class BulkProcess:
    """bulk process for datafilelist based on indexfilelist

    Attributes:
        overall_system  :   final good system
        overall_sysbad  :   all bad system

        overall_prob_initial    :   based on read data file
        overall_prob_final      :   based on read data file

    Note:
        fragments is input human-readable number, starting at 1
    """
    def __init__ (self,datafilelist,*args,**kwargs):
        fd = []
        for f in datafilelist:
            if os.path.isfile(f):
                fd.append(f)
            else:
                print('Warning: not a data file < {:} >, ignoring'.format(f))
        self.datafilelist = fd

        self.indexfilelist = []
        if 'indexfilelist' in kwargs and kwargs['indexfilelist'] is not None:
            fd = []
            for f in kwargs['indexfilelist']:
                if os.path.isfile(f):
                    fd.append(f)
                else:
                    print('Warning: not an index file < {:} >, ignoring'.format(f))
            self.indexfilelist = fd
        
        self.bool_save_images = True
        if 'bool_save_images' in kwargs:
            self.bool_save_images = kwargs['bool_save_images']
        self.bool_save_images = True if self.bool_save_images is True else False

        if len(self.datafilelist) == 0:
            print('Warning: no inputs')
            raise ValueError('no inputs')

        self.args = args
        self.kwargs = kwargs

        self.overall_system = []
        self.overall_sysbad = []
        self.overall_prob_initial = []
        self.overall_prob_final = []

        self.overall_energy = []

        self.force_double_check = False
        if 'force_double_check' in kwargs:
            if kwargs['force_double_check'] is True:
               self.force_double_check = True



    def run(self):
        print('Note: reading all data file inputs')
        systemlist,energylist = self.get_datalist(self.datafilelist)
        self.molnms = [len(i) for i in systemlist]

        sysndxlist = []
        if len(self.indexfilelist) != 0:
            print('Note: reading all index file inputs')
            sysndxlist,tmp = self.get_datalist(self.indexfilelist)

        # connections only need to be calculated once
        self.get_connections(system=systemlist[0][0])

        # prompt for double check
        bo = True
        if self.force_double_check:
            print('\nCheck: data files:')
            for cnt,fd in enumerate(self.datafilelist):
                print('   => {:} -- molnms {:}'.format(fd,self.molnms[cnt]))

            if len(self.indexfilelist) != 0:
                print('Check: index filesd:')
                for fd in self.indexfilelist:
                    print('   => {:}'.format(fd))

            ltmp = self.kwargs['fragments']
            if ltmp is not None and len(ltmp) != 0:
                print('Check: molecule fragments:')
                lt = []
                for i in ltmp: lt.append([j+1 for j in i])
                print('   => {:}'.format(lt))

            ltmp = self.kwargs['bcon']
            if ltmp is not None and len(ltmp) != 0:
                print('Check: bond connection:')
                lt = []
                for i in ltmp: lt.append([j+1 for j in i])
                print('   => {:}'.format(lt))

            ltmp = self.kwargs['acon']
            if ltmp is not None and len(ltmp) != 0:
                print('Check: angle connection:')
                lt = []
                for i in ltmp: lt.append([j+1 for j in i])
                print('   => {:}'.format(lt))
            
            print('Check: total inputs < {:} >'.format(sum(self.molnms)))
            
            stmp = 'ON' if self.kwargs['obonds'] else 'OFF'
            print('Check: bonds probability is < {:} >'.format(stmp))

            stmp = 'ON' if self.kwargs['oangles'] else 'OFF'
            print('Check: angles probability is < {:} >'.format(stmp))

            stmp = 'ON' if self.kwargs['opar'] else 'OFF'
            print('Check: par probability < {:} > (time consuming)'.format(stmp))

            stmp = 'ON' if self.kwargs['oall'] else 'OFF'
            print('Check: all probability is < {:} >'.format(stmp))

            stmp = 'ON' if self.kwargs['bool_save_images'] else 'OFF'
            print('Check: image outputs are < {:} >'.format(stmp))

            print('\nDo you want to continue? y/yes, else not. Input: ',end='')
            tmp = input().lower()
            if tmp not in ['y','yes']: bo = False
            if not bo:
                print('Note: you decided to quit, nothing will be processed')
            print()

        if bo:
            self.filtration_on_self(systemlist,energylist)
            self.filtration_on_cross(self.overall_system, sysndxlist)
            if len(self.datafilelist) > 1: self.filtration_on_overall()
            self.file_print()
            self.save_data()



    def file_print(self):
        print('\nNote: saving bulk process results ...')
        outmolnms = len(self.overall_system)
        print('Note: final molnms: < {:} >'.format(outmolnms))
        outratio = round(outmolnms/sum(self.molnms),2)
        print('Note: filtration ratio: < {:} >'.format(outratio))

        if 'fname' in self.kwargs and self.kwargs['fname'] is not None:
            if len(self.kwargs['fname'].split()) != 0:
                fd = self.kwargs['fname']
            else:
                fd = 'system.txt'
        else:
            fd = 'system.txt'
        
        ext = 'txt'if self.kwargs['ftype'] is None else self.kwargs['ftype']
        self.kwargs['fname'] = file_gen_new(fd,fextend=ext)

        self.kwargs['energy'] = self.overall_energy
        FD = SaveFile(self.overall_system,*self.args,**self.kwargs)
        FD.run()

        if self.bool_save_images:
            fimgs = {'bonds':[], 'angles':[]}
            for cnt,initial in enumerate(self.overall_prob_initial):
                # save images for bonds
                fout = {'par':[], 'all':None}
                if len(initial['bonds']) != 0:
                    fg = 'bulk-image-bonds'
                    key = 'bonds'
                    ini = initial['bonds']
                    fin = self.overall_prob_final[cnt]['bonds']
                    dt = self.tolerance[0]
                    fout = self.save_images(ini,fin,fg,dt=dt,key=key)
                fimgs['bonds'].append(fout)

                # save images for angles
                fout = {'par':[], 'all':None}
                if len(initial['angles']) != 0:
                    fg = 'bulk-image-angles'
                    key = 'angles'
                    ini = initial['angles']
                    fin = self.overall_prob_final[cnt]['angles']
                    dt = self.tolerance[1]
                    fout = self.save_images(ini,fin,fg,dt=dt,key=key)
                fimgs['angles'].append(fout)
        
        ftot = file_gen_new('bulk-process-info')
        print('\nNote: check summary file for more info: < {:} >'.format(ftot))
        with open(ftot,'wt') as f:
            f.write('Note: bulk process for input files:\n')
            for cnt,fd in enumerate(self.datafilelist):
                f.write('  => {:} -- molnms {:}\n'.format(fd,self.molnms[cnt]))
            f.write('\nNote: index files:\n')
            for fd in self.indexfilelist:
                f.write('  => {:}\n'.format(fd))
            f.write('\nNote: result file:\n')
            f.write('  => {:} -- molnms {:}\n'.format(FD.fname,outmolnms))
            f.write('\nNote: filtration ratio: {:}\n'.format(outratio))

            bcon = self.kwargs['bcon']
            acon = self.kwargs['acon']
            f.write('\nNote: bonds connections:\n')
            f.write('  => {:}\n'.format(bcon))
            f.write('\nNote: angles connections:\n')
            f.write('  => {:}\n'.format(acon))

            f.write('\nNote: tolerance:\n')
            f.write('  => {:}\n'.format(self.tolerance))

            if not self.bool_save_images: return

            for ndx,df in enumerate(self.datafilelist):
                f.write('\n\nNote: images for file: < {:} >\n'.format(df))
                imgs = fimgs['bonds'][ndx]
                if len(imgs) != 0:
                    if imgs['all'] is not None:
                        f.write('Note: bonds images for all in pot of stew:\n')
                        f.write('  => all --> {:}\n'.format(imgs['all']))
                    if len(imgs['par']) != 0:
                        f.write('Note: bonds images for each par entry:\n')
                        for cnt,fd in enumerate(imgs['par']):
                            f.write('  => {:} --> {:}\n'.format(bcon[cnt],fd))

                imgs = fimgs['angles'][ndx]
                if len(imgs) != 0:
                    if imgs['all'] is not None:
                        f.write('Note: angles images for all in pot of stew:\n')
                        f.write('  => all --> {:}\n'.format(imgs['all']))
                    if len(imgs['par']) != 0:
                        f.write('Note: angles images for each par entry:\n')
                        for cnt,fd in enumerate(imgs['par']):
                            f.write('  => {:} --> {:}\n'.format(acon[cnt],fd))

            # if number of data file inpust bigger than 1, overall filtration
            if len(self.datafilelist) > 1:
                f.write('\n\nNote: overall filtration on all files\n')
                imgs = fimgs['bonds'][-1]
                if len(imgs) != 0:
                    if imgs['all'] is not None:
                        f.write('Note: bonds images for overall in pot of stew:\n')
                        f.write('  => overall --> {:}\n'.format(imgs['all']))
                    if len(imgs['par']) != 0:
                        f.write('Note: bonds images for each par entry:\n')
                        for cnt,fd in enumerate(imgs['par']):
                            f.write('  => {:} --> {:}\n'.format(bcon[cnt],fd))

                imgs = fimgs['angles'][-1]
                if len(imgs) != 0:
                    if imgs['all'] is not None:
                        f.write('Note: angles images for overall in pot of stew:\n')
                        f.write('  => overall --> {:}\n'.format(imgs['all']))
                    if len(imgs['par']) != 0:
                        f.write('Note: angles images for each par entry:\n')
                        for cnt,fd in enumerate(imgs['par']):
                            f.write('  => {:} --> {:}\n'.format(acon[cnt],fd))



    def save_data(self):
        def gen_outputs(ini,fin,ctype,key):
            out = ''
            itxt = '@{:}-INI-{:}'.format(ctype.upper(),key.upper())
            ftxt = '@{:}-FIN-{:}'.format(ctype.upper(),key.upper())
            for m,p in enumerate(ini):
                if len(p) == 0:
                    out += '{:}-X-{:}   NONE\n'.format(itxt,m+1)
                    out += '{:}-Y-{:}   NONE\n'.format(itxt,m+1)
                else:
                    out += '{:}-X-{:}   {:}\n'.format(itxt,m+1,p[1])
                    out += '{:}-Y-{:}   '.format(itxt,m+1)
                    for s in p[0]: out += '{:}   '.format(s)
                    out += '\n'
                u = fin[m]
                if len(u) == 0:
                    out += '{:}-X-{:}   NONE\n'.format(ftxt,m+1)
                    out += '{:}-Y-{:}   NONE\n'.format(ftxt,m+1)
                else:
                    out += '{:}-X-{:}   {:}\n'.format(ftxt,m+1,u[1])
                    out += '{:}-Y-{:}   '.format(ftxt,m+1)
                    for s in u[0]: out += '{:}   '.format(s)
                    out += '\n'
            return out


        fdata = file_gen_new('bulk-probability-data')
        print('Note: probability data is saved to < {:} >'.format(fdata))
        with open(fdata,'wt') as f:
            f.write('@TOLERANCE   {:}   {:}\n\n\n'.format(*self.tolerance))
            for cnt,initial in enumerate(self.overall_prob_initial):
                if cnt < len(self.datafilelist):
                    finfo = self.datafilelist[cnt] 
                else:
                    finfo = 'OVERALL'
                f.write('@FILE {:}\n\n'.format(finfo))

                # for bonds
                if len(initial['bonds']) == 0:
                    f.write('@BONDS   NONE\n\n')
                elif len(initial['bonds']) != 0:
                    ini = initial['bonds']
                    fin = self.overall_prob_final[cnt]['bonds']

                    ipar = ini[0]
                    fpar = fin[0]
                    if len(ipar) == 0:
                        f.write('@BONDS-PAR   NONE\n')
                    else:
                        out = gen_outputs(ipar,fpar,'BONDS','PAR')
                        f.write(out)
                    f.write('\n')

                    new = ''
                    iall = ini[1]
                    if len(iall) == 0:
                        f.write('@BONDS-INI-ALL-X   NONE\n')
                        f.write('@BONDS-INI-ALL-Y   NONE\n')
                    else:
                        new += '@BONDS-INI-ALL-X   {:}\n'.format(iall[1])
                        new += '@BONDS-INI-ALL-Y   '
                        for s in iall[0]: new += '{:}   '.format(s)
                        new += '\n'
                    fall = fin[1]
                    if len(fall) == 0:
                        f.write('@BONDS-FIN-ALL-X   NONE\n')
                        f.write('@BONDS-FIN-ALL-Y   NONE\n')
                    else:
                        new += '@BONDS-FIN-ALL-X   {:}\n'.format(fall[1])
                        new += '@BONDS-FIN-ALL-Y   '
                        for s in fall[0]: new += '{:}   '.format(s)
                        new += '\n'
                    f.write(new)
                    f.write('\n')
                
                # for angles
                if len(initial['angles']) == 0:
                    f.write('@ANGLES   NONE\n\n')
                elif len(initial['angles']) != 0:
                    ini = initial['angles']
                    fin = self.overall_prob_final[cnt]['angles']

                    ipar = ini[0]
                    fpar = fin[0]
                    if len(ipar) == 0:
                        f.write('@ANGLES-PAR   NONE\n')
                    else:
                        out = gen_outputs(ipar,fpar,'ANGLES','PAR')
                        f.write(out)
                    f.write('\n')

                    new = ''
                    iall = ini[1]
                    if len(iall) == 0:
                        f.write('@ANGLES-INI-ALL-X   NONE\n')
                        f.write('@ANGLES-INI-ALL-Y   NONE\n')
                    else:
                        new += '@ANGLES-INI-ALL-X   {:}\n'.format(iall[1])
                        new += '@ANGLES-INI-ALL-Y   '
                        for s in iall[0]: new += '{:}   '.format(s)
                        new += '\n'
                    fall = fin[1]
                    if len(fall) == 0:
                        f.write('@ANGLES-FIN-ALL-X   NONE\n')
                        f.write('@ANGLES-FIN-ALL-Y   NONE\n')
                    else:
                        new += '@ANGLES-FIN-ALL-X   {:}\n'.format(fall[1])
                        new += '@ANGLES-FIN-ALL-Y   '
                        for s in fall[0]: new += '{:}   '.format(s)
                        new += '\n'
                    f.write(new)
                    f.write('\n')
                f.write('\n\n')



    def save_images(self,ini,fin,fg,dt=None,key=None):
        """separately save images

        Return:
            fout:   generated files name dict

            key :   None means file not exist
                par :   1D  :   List[str]
                all :   str
        """
        fout = {'par':[], 'all':None}
        # for par
        if len(ini[0]) != 0:
            ftmp = fg + '-par'
            for ndx,par in enumerate(ini[0]):
                fgp = file_gen_new(ftmp,fextend='jpg',foriginal=False)
                fbo = plot_filtration_save_image(
                    par,
                    fin[0][ndx],
                    dt=dt,
                    fname=fgp,
                    key=key,
                )
                fout['par'].append(fgp if fbo else None)
        # for all
        if len(ini[1]) != 0:
            fga = fg + '-all'
            fga = file_gen_new(fga,fextend='jpg',foriginal=False)
            fbo = plot_filtration_save_image(
                ini[1],
                fin[1],
                dt=dt,
                fname=fga,
                key=key,
            )
            fout['all'] = fga if fbo else None
        return fout



    def get_datalist(self,filelist):
        datalist = []
        energylist = []
        for f in filelist:
            DATA = ReadFile(f)
            DATA.run()
            print('Note: for file < {:} >, number of inputs < {:} >'.format(f,len(DATA.system)))
            datalist.append(DATA.system)
            energylist.append(DATA.energy)
        return datalist,energylist



    def get_connections(self,system):
        """
        Rule:
            There are three type of inputs,

            1) bcon         ->  bob True if has else False
            2) acon         ->  boa True if has else False
            3) fragments    ->  bog True if has else False

            now, we think as user input is always on the first priority,
            and bcon & acon at the second, fragments in the last

            thus fragments is an optional input, which only works;
            when itself as the input and when bcon or acon not exist

            a special case needs to be considered, when all three does not
            exist, then the connections are calculated by the input system

            after process, all three attributes will always be attached
        """
        bob = False
        if 'bcon' in self.kwargs and self.kwargs['bcon'] is not None:
            bob = True
        boa = False
        if 'acon' in self.kwargs and self.kwargs['acon'] is not None:
            boa = True
        bog = False
        if 'fragments' in self.kwargs and self.kwargs['fragments'] is not None:
            bog = True
        if not bog:
            # to simplify, make fragments always exist
            BOND = BondPerception(system,3)
            BOND.run()
            self.kwargs['fragments'] = BOND.fragments

        # only check user input fragments
        if bog:
            ux = max([max(i) for i in self.kwargs['fragments']])
            if ux > len(system):
                print('Warning: atom index is too big: < {:} >'.format(ux))
                raise ValueError('wrongly defined')
            self.kwargs['fragments'] = self.check_user_input_connections(
                self.kwargs['fragments'],
                self.kwargs['fragments'],
            )

        bcon,acon = func_calc_connection(self.kwargs['fragments'])
        # take care of special case, when input is None
        if bob:
            self.kwargs['bcon'] = self.check_user_input_connections(
                self.kwargs['bcon'],
                self.kwargs['fragments'],
            )
        else:
            if 'bcon' in self.kwargs and self.kwargs['bcon'] is None:
                self.kwargs['bcon'] = bcon
            else:
                self.kwargs['bcon'] = []
        if boa:
            self.kwargs['acon'] = self.check_user_input_connections(
                self.kwargs['acon'],
                self.kwargs['fragments'],
            )
        else:
            if 'acon' in self.kwargs and self.kwargs['acon'] is None:
                self.kwargs['acon'] = acon
            else:
                self.kwargs['acon'] = []



    def check_user_input_connections(self,ul,fl):
        for ndx in ul:
            if len(set(ndx)) != len(ndx):
                print('Warning: repeats at: < {:} >'.format(ndx))
                raise ValueError('wrongly defined')

        ux = max([max(i) for i in ul])
        fx = max([max(i) for i in fl])
        if ux-1 > fx:
            print('Warning: atom index exceeding: < {:} >'.format(ux))
            raise ValueError('wrongly defined')
        # convert to python readable number, starts at 0
        ptmp = []
        for i in ul: ptmp.append([j-1 for j in i])
        for i in ptmp:
            if min(i) < 0:
                print('Warning: atom index starts at 1')
                raise ValueError('wrongly defined')
        return ptmp



    def filtration_on_cross(self,system,sysndxlist):
        for cnt,sysndx in enumerate(sysndxlist):
            print('Note: cross filtration on index file < {:} >'.format(self.indexfilelist[cnt]))
            CFF = CrossFiltration(system, sysndx=sysndx, *self.args, **self.kwargs)

            self.overall_system = CFF.system
            for mol in CFF.sysbad: self.overall_sysbad.append(mol)

            self.overall_energy = [e for i,e in enumerate(self.overall_energy) if i not in CFF.reflist]



    def filtration_on_self(self,datalist,energylist):
        for cnt,system in enumerate(datalist):
            print('Note: self filtration for file < {:} >'.format(self.datafilelist[cnt]))
            FF = Filtration(system, *self.args, **self.kwargs)
            FF.run()
            self.tolerance = FF.tolerance

            for mol in FF.system: self.overall_system.append(mol)
            for mol in FF.sysbad: self.overall_sysbad.append(mol)

            self.overall_prob_initial.append(FF.prob_initial)
            self.overall_prob_final.append(FF.prob_final)

            # update energy
            energy = energylist[cnt]
            for i,e in enumerate(energy):
                if i not in FF.reflist:
                    self.overall_energy.append(e)



    def filtration_on_overall(self):
        print('Note: overall filtration on all filtered inputs')
        FF = Filtration(self.overall_system, *self.args, **self.kwargs)
        FF.run()

        self.overall_system = FF.system
        for mol in FF.sysbad: self.overall_sysbad.append(mol)

        self.overall_prob_initial.append(FF.prob_initial)
        self.overall_prob_final.append(FF.prob_final)

        self.overall_energy = [e for i,e in enumerate(self.overall_energy) if i not in FF.reflist]




def parsecmd():
    """Parse command line input"""
    def parse_remove_chars(line):
        line = line.replace('[',' ').replace(']',' ').replace(';',',')
        return line
    
    def parse_in_line(line,bobcon=False,boacon=False):
        line = parse_remove_chars(line)
        lt = line.split(',')
        ls = []
        for t in lt:
            lp = t.split()
            if len(lp) == 0: continue
            if bobcon:
                if len(lp) == 2:
                    ls.append(lp)
                else:
                    print('Warning: wrong bond connection < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            elif boacon:
                if len(lp) == 3:
                    ls.append(lp)
                else:
                    print('Warning: wrong angle connection < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            else:
                ls.append(lp)

        if len(ls) == 0: return []
        
        lm = []
        for t in ls:
            lx = []
            for v in t:
                try:
                    lx.append(int(v))
                except ValueError:
                    print('Warning: wrong input < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            lm.append(lx)
        
        return lm


    parser = argparse.ArgumentParser(
        description='Conformation Filtration',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION,
    )
    parser.add_argument(
        '-f','--datafiles',
        help='Data files, separate by space or comma',
        nargs='+',
    )
    parser.add_argument(
        '-d','--indexfiles',
        help='Index files, as the filtration reference',
        nargs='+',
    )
    parser.add_argument(
        '-t','--tolerance',
        help='In sequence, Bonds(Angstrom), Angles(Degree)',
        nargs='+',
        metavar='inc',
    )
    parser.add_argument(
        '-b','-bcon','--bcon',
        help='Bond Connections, in pairs, separate by comma',
        nargs='+',
        metavar='B',
    )
    parser.add_argument(
        '-a','-acon','--acon',
        help='Angle Connections, in pairs, separate by comma',
        nargs='+',
        metavar='A',
    )
    parser.add_argument(
        '-g','--fragments',
        help='Molecular fragments, separate by comma',
        nargs='+',
        metavar='G',
    )
    parser.add_argument(
        '-ob','-obonds','--obonds',
        help='Turn on bonds probability calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-oa','-oangles','--oangles',
        help='Turn on angles probability calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-op','-opar','--opar',
        help='Turn on calculation for each par, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '--no-oall',
        help='Turn off calculation for in all, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-nc','--no-force-double-check',
        help='Turn off double check before execution',
        action='store_true',
    )
    parser.add_argument(
        '-ni','--no-image-outputs',
        help='Turn off image outputs',
        action='store_true',
    )
    parser.add_argument(
        '-at','--force-real-atom-type',
        help='Outputs using real atom type, if possible',
        action='store_true',
    )
    parser.add_argument(
        '--features',
        help='Show develop features',
        action='store_true',
    )
    parser.add_argument(
        '-p','--file-format-explanations',
        help='Show input file format explanations',
        action='store_true',
    )
    parser.add_argument(
        '-o','--fname',
        help='Output file name',
    )
    parser.add_argument(
        '-ft','--ftype',
        help='Output file type, default txt',
    )

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    args = parser.parse_args(sys.argv[1:])
    if args.features:
        for i in FEATURES:
            print(i)
        exit()
    
    if args.file_format_explanations:
        print(FILEFORMAT)
        exit()

    # default settings
    fdict = {
        #'datafilelist'          :   None,
        'indexfilelist'         :   None,
        'bcon'                  :   None,
        'acon'                  :   None,
        'fragments'             :   None,
        'force_double_check'    :   True,
        'force_real_atom_type'  :   None,
        'bool_save_images'      :   True,
        'tolerance'             :   None,
        'fname'                 :   None,
        'ftype'                 :   'txt',
        'obonds'                :   None,
        'oangles'               :   None,
        'opar'                  :   None,
        'oall'                  :   True,
    }
    if args.datafiles is None:
        print('Warning: -f/--datafiles is missing')
        exit()

    stmp = parse_remove_chars(' '.join(args.datafiles)).replace(',',' ')
    datafilelist = stmp.split()

    if args.indexfiles is not None:
        stmp = parse_remove_chars(' '.join(args.indexfiles)).replace(',',' ')
        fdict['indexfilelist'] = stmp.split()

    if args.bcon is not None:
        fdict['bcon'] = parse_in_line(' '.join(args.bcon),bobcon=True)
    if args.acon is not None:
        fdict['acon'] = parse_in_line(' '.join(args.acon),boacon=True)
    if args.fragments is not None:
        fdict['fragments'] = parse_in_line(' '.join(args.fragments))

    if args.obonds: fdict['obonds'] = True
    if args.oangles: fdict['oangles'] = True
    if args.opar: fdict['opar'] = True
    if args.no_oall: fdict['oall'] = False
    if args.no_force_double_check: fdict['force_double_check'] = False
    if args.force_real_atom_type: fdict['force_real_atom_type'] = True
    if args.no_image_outputs: fdict['bool_save_images'] = False

    if args.fname is not None: fdict['fname'] = args.fname
    if args.ftype is not None: fdict['ftype'] = args.ftype

    if args.tolerance is not None:
        stmp = parse_remove_chars(' '.join(args.tolerance)).replace(',',' ')
        lt = stmp.split()
        bo = False
        if len(lt) == 0:
            pass
        elif len(lt) == 2:
            try:
                t1 = float(lt[0])
                t2 = float(lt[1])
                fdict['tolerance'] = [t1,t2]
            except ValueError:
                bo = True
        else:
            bo = True
        if bo:
            print('Warning: wrong defined < {:} >'.format(args.tolerance))
            raise ValueError('wrongly defined')

    FF = BulkProcess(datafilelist,**fdict,)
    FF.run()




if __name__ == '__main__':
    parsecmd()

