#!/usr/bin/env python3

import io
import os
import numpy as np
import itertools
import matplotlib.pyplot as plt

# for testing
import tempfile

FEATURES = [
    'version 0.1    : MOlecular Structure Sampling',
    'version 0.2    : add ClosedForm',
    'version 0.3    : add unified ParConfig',
    'version 0.4    : add func translation',
]


VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


FILEFORMAT = """
Input File Format for BOSS Output Filtration

Currently, only following file formats are supported.

For txt file:

    molecules are separated by new line, line starts with char '#' will
    be ignored, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   # x x x [=] energy
            |   atom-1   x    y    z
    mol     |   atom-2   x    y    z
            |   atom-3   x    y    z
            |   ...
            [   <new line>


For xsf file:

    molecules are separated by '#', keyword 'ATOMS' is important,
    case sensitive, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   # x x x [=] energy
            |
            |   ATOMS
    mol     |   atom-1   x    y    z    else
            |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...


For xyz file:

    molecules are separated by new line, line starts with char '#' will
    be ignored, if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the first line,
          equal sign, =, can be used.

            [   #  x x x [=] energy
            |   atom-1   x    y    z    else
    mol     |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...


For com file:

    molecules are separated by '#'
    if any errors happen, the whole set will be skipped.

    note: energy info has to be at the end of the title line,
          equal sign, =, can be used.

            [   %cmd=..
            |   %cmd=...
            |   %cmd=....
            |   # basis_set
            |
    mol     |   title   x x x [=] energy
            |
            |   charge_spin
            |   atom-1   x    y    z    else
            |   atom-2   x    y    z    else
            |   atom-3   x    y    z    else
            [   ...

"""



class ReadFileMultiple:
    """
    Initialization:
        file :
            format:
                txt file
                xsf file
                xyz file

    Attributes:
        system : System[Mol[Atom, ...]] : List[List[[atomtype, x,y,z], ...]]
    """
    def __init__(self,file,ext=None):
        self.nice = True
        self.info = ''
        self.file = file

        # decide file format
        if ext is None:
            if isinstance(file,int): file = str(file)
            ndx = file.rfind('.')
            if ndx == -1 or ndx + 1 >= len(file):
                self.ext = None
            else:
                self.ext = file[ndx+1:].lower()
            if self.ext is None: self.ext = 'txt'
        else:
            self.ext = ext

        if self.ext not in ['txt','xsf','xyz']:
            self.nice = False
            self.info = 'Warning: not support: file format: {:}'.format(ext)
            print(self.info)
            return



    def run(self):
        print('Note: reading file < {:} > ...'.format(self.file))
        fn = getattr(self,'read_'+self.ext)
        fn()



    def read_xsf(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        promol = []
        i = 0
        while i < len(profile):
            line = profile[i]
            if line.find('#') != -1:
                j = i + 1
                ls = [[line.strip(),i], ]
                while j < len(profile):
                    sub = profile[j]
                    if sub.find('#') != -1: break
                    sub = sub.strip()
                    if len(sub) != 0: ls.append([sub,j])
                    j += 1
                promol.append(ls)
                i = j
            else:
                i += 1

        self.system = []
        for cnt,mol in enumerate(promol):
            bo = False
            if len(mol) <= 2:
                bo = True
                errnum = mol[0][1] + 1
                errline = 'Wrong format'
            else:
                if mol[0][0][0] != '#' == -1 or mol[1][0] != 'ATOMS':
                    bo = True
                    if mol[0][0][0] != '#':
                        errnum = mol[0][1] + 1
                    else:
                        errnum = mol[1][1] + 1
                    errline = 'Wrong format'

            if not bo:
                ls = []
                for t in mol[2:]:
                    # atom info
                    atom = t[0].split()
                    if len(atom) >= 4:
                        try:
                            x = float(atom[1])
                            y = float(atom[2])
                            z = float(atom[3])
                        except ValueError:
                            bo = True
                            errnum = t[1] + 1
                            errline = t[0]
                    else:
                        bo = True
                        errnum = t[1] + 1
                        errline = t[0]
                    if bo:
                        break
                    else:
                        ls.append([atom[0],x,y,z])

            if bo:
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                self.system.append(ls)



    def read_txt(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        # List[List[[atomtype, x, y, z], ...]]
        promol = []
        mol = []
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

        self.system = []
        for mol in promol:
            if mol[0][0][0] == '#': mol = mol[1:]

            bo = False
            ls = []
            for t in mol:
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
            if bo:
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                self.system.append(ls)



    def read_xyz(self):
        with open(self.file,mode='rt') as f:
            profile = f.readlines()

        # List[List[[atomtype, x, y, z], ...]]
        promol = []
        mol = []
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

        self.system = []
        for mol in promol:
            bo = False
            if len(mol) <= 2: bo = True
            if not bo:
                try:
                    atomnum = int(mol[0][0])
                    if atomnum + 2 != len(mol):
                        raise ValueError
                except ValueError:
                    bo = True

            if bo:
                errnum = mol[0][1] + 1
                errline = mol[0][0]
            else:
                ls = []
                for t in mol[2:]:
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

            if bo:
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                self.system.append(ls)




def test_class_ReadFileMultiple():
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000
        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325

        H           -1.081          -1.003          1.416
        H            0.576           -0.336         1.933
        C           -0.329            1.632         2.341
        O            1.268            0.242         1.064
        N           -0.017            0.349         4.373
        10           1.331            2.040         3.305
        H            0.155           -0.527         2.467
        H            1.118            0.107         3.318
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()

    for mol in RF.system:
        for at in mol:
            print(at)
        print()
        print()
    

    fxsf = """
        #
        ATOMS

        C         -0.22451        -0.34646         1.16745
        C         -0.91697         0.72807         1.99904
        O              0.0             0.0             0.0
        N         -0.97796        -0.00661         3.14107
        O         -1.95042        -0.78556         3.32521
        O             -0.0            -0.0         3.93472


        #
        ATOMS

        C         -0.23085        -0.33801         1.17685
        C         -0.91195          0.7211         1.99194
        O         -0.00509         0.00579         0.00885
        N         -0.97441        -0.01109         3.13549
        O         -1.94755        -0.78899         3.32046
        O          0.00286        -0.00343         3.92997
        H         -1.83972         0.81292         1.44776
        H         -0.15568          1.4912         1.97403
        H         -1.08767        -0.99342         1.42623
        H          0.56856        -0.32597         1.94281
    """

    print('\n\n')
    fp = None
    fp = tempfile.TemporaryFile()
    fp.write(fxsf.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name,ext='xsf')
    # fp will be destoried inside open!
    RF.run()

    for mol in RF.system:
        for at in mol:
            print(at)
        print()
        print()


    fxyz = """
        10
        Properties=species:S:1:pos:R:3 energy=-223995.61806643562
        C      -0.16801     -0.53291      1.19682
        C       -0.6905      0.48737      2.20263
        O          -0.0          0.0          0.0
        N      -0.89651     -0.06647      3.58017
        O      -1.97078     -0.55868       3.8878
        O          -0.0         -0.0      4.40678
        H      -1.63684      0.92363      1.79682
        H       0.01307      1.35677      2.21605
        H      -0.89021      -1.4027      1.19644
        H       0.78356     -0.96844      1.62482
    
    """

    print('\n\n')
    fp = None
    fp = tempfile.TemporaryFile()
    fp.write(fxyz.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name,ext='xyz')
    # fp will be destoried inside open!
    RF.run()

    for mol in RF.system:
        for at in mol:
            print(at)
        print()
        print()




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




class AtomInfo:
    """
    Method:
        get_atom    : return Atom(s,n,r,m,name,xyz)

    Args:
        s (str)     : sign, atomtype either can be real atom type or its number
        n (int)     : number in periodic table
        r (float)   : radius    (Angstrom)
        m (float)   : mass  (g/mol)
        name (str)  : full name
        xyz (List)  : cartesian coordinates (Angstrom), default [0., 0., 0.]

    Reference:
        Visualizing atomic sizes and molecular shapes with the classical
        turning surface of the Kohn–Sham potential

        Egor Ospadov, Jianmin Tao, Viktor N. Staroverov, John P. Perdew

        Proceedings of the National Academy of Sciences Dec 2018,
        115 (50) E11578-E11585; DOI: 10.1073/pnas.1814300115


        Atomic radius got from:

        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Note:
        A special dummy atom X is defined:
            sign=X, number=0, radius=0.0, mass=1.0, name=Dummy

        If conflict happens, Atom is defined in the sequence:
            1) sign
            2) number
            3) name

        If no inputs, dummy atom will be returned

        If any errors happen, dummy atom will be returned
    """
    ATOM_PROPERTY = """sign,number,radius,mass,name
        X,0,0,1,Dummy
        H,1,0.23,1.00794,Hydrogen
        He,2,0.93,4.002602,Helium
        Li,3,0.68,6.941,Lithium
        Be,4,0.35,9.01218,Beryllium
        B,5,0.83,10.811,Boron
        C,6,0.68,12.011,Carbon
        N,7,0.68,14.00674,Nitrogen
        O,8,0.68,15.9994,Oxygen
        F,9,0.64,18.998403,Fluorine
        Ne,10,1.12,20.1797,Neon
        Na,11,0.97,22.989768,Sodium
        Mg,12,1.1,24.305,Magnesium
        Al,13,1.35,26.981539,Aluminum
        Si,14,1.2,28.0855,Silicon
        P,15,1.05,30.973762,Phosphorus
        S,16,1.02,32.066,Sulfur
        Cl,17,0.99,35.4527,Chlorine
        Ar,18,1.57,39.948,Argon
        K,19,1.33,39.0983,Potassium
        Ca,20,0.99,40.078,Calcium
        Sc,21,1.44,44.95591,Scandium
        Ti,22,1.47,47.88,Titanium
        V,23,1.33,50.9415,Vanadium
        Cr,24,1.35,51.9961,Chromium
        Mn,25,1.35,54.93805,Manganese
        Fe,26,1.34,55.847,Iron
        Co,27,1.33,58.9332,Cobalt
        Ni,28,1.5,58.6934,Nickel
        Cu,29,1.52,63.546,Copper
        Zn,30,1.45,65.39,Zinc
        Ga,31,1.22,69.723,Gallium
        Ge,32,1.17,72.61,Germanium
        As,33,1.21,74.92159,Arsenic
        Se,34,1.22,78.96,Selenium
        Br,35,1.21,79.904,Bromine
        Kr,36,1.91,83.8,Krypton
        Rb,37,1.47,85.4678,Rubidium
        Sr,38,1.12,87.62,Strontium
        Y,39,1.78,88.90585,Yttrium
        Zr,40,1.56,91.224,Zirconium
        Nb,41,1.48,92.90638,Niobium
        Mo,42,1.47,95.94,Molybdenum
        Tc,43,1.35,97.9072,Technetium
        Ru,44,1.4,101.07,Ruthenium
        Rh,45,1.45,102.9055,Rhodium
        Pd,46,1.5,106.42,Palladium
        Ag,47,1.59,107.8682,Silver
        Cd,48,1.69,112.411,Cadmium
        In,49,1.63,114.818,Indium
        Sn,50,1.46,118.71,Tin
        Sb,51,,121.76,Antimony
        Te,52,1.47,127.6,Tellurium
        I,53,1.4,126.90447,Iodine
        Xe,54,1.98,131.29,Xenon
        Cs,55,1.67,132.90543,Cesium
        Ba,56,1.34,137.327,Barium
        La,57,1.87,138.9055,Lanthanum
        Ce,58,1.83,140.115,Cerium
        Pr,59,1.82,140.90765,Praseodymium
        Nd,60,1.81,144.24,Neodymium
        Pm,61,1.8,144.9127,Promethium
        Sm,62,1.8,150.36,Samarium
        Eu,63,1.99,151.965,Europium
        Gd,64,1.79,157.25,Gadolinium
        Tb,65,1.76,158.92534,Terbium
        Dy,66,1.75,162.5,Dysprosium
        Ho,67,1.74,164.93032,Holmium
        Er,68,1.73,167.26,Erbium
        Tm,69,1.72,168.93421,Thulium
        Yb,70,,173.04,Ytterbium
        Lu,71,,174.967,Lutetium
        Hf,72,,178.49,Hafnium
        Ta,73,,180.9479,Tantalum
        W,74,,183.84,Tungsten
        Re,75,,186.207,Rhenium
        Os,76,,190.23,Osmium
        Ir,77,,192.22,Iridium
        Pt,78,,195.08,Platinum
        Au,79,,196.96654,Gold
        Hg,80,,200.59,Mercury
        Tl,81,,204.3833,Thallium
        Pb,82,,207.2,Lead
        Bi,83,,208.98037,Bismuth
        Po,84,,208.9824,Polonium
        At,85,,209.9871,Astatine
        Rn,86,,222.0176,Radon
        Fr,87,,223.0197,Francium
        Ra,88,,226.0254,Radium
        Ac,89,,227.0278,Actinium
        Th,90,,232.0381,Thorium
        Pa,91,,231.03588,Protactinium
        U,92,,238.0289,Uranium
        Np,93,,237.048,Neptunium
        Pu,94,,244.0642,Plutonium
        Am,95,,243.0614,Americium
        Cm,96,,247.0703,Curium
        Bk,97,,247.0703,Berkelium
        Cf,98,,251.0796,Californium
        Es,99,,252.083,Einsteinium
        Fm,100,,257.0951,Fermium
        Md,101,,258.1,Mendelevium
        No,102,,259.1009,Nobelium
        Lr,103,,262.11,Lawrencium
        Rf,104,,261,Rutherfordium
        Db,105,,262,Dubnium
        Sg,106,,266,Seaborgium
        Bh,107,,264,Bohrium
        Hs,108,,269,Hassium
        Mt,109,,268,Meitnerium
        Ds,110,,269,Darmstadtium
        Rg,111,,272,Roentgenium
        Cn,112,,277,Copernicium
        Uut,113,,,Ununtrium
        Fl,114,,289,Flerovium
        Uup,115,,,Ununpentium
        Lv,116,,,Livermorium
        Uus,117,,,Ununseptium
        Uuo,118,,,Ununoctium"""


    def __init__(self):
        self.atominfo = self.pro_atom_property(self.ATOM_PROPERTY)



    class Atom:
        def __init__(self,s,n,r,m,name,xyz):
            self.s = s
            self.n = n
            self.r = r
            self.m = m
            self.name = name
            self.xyz = xyz



    def copy_atom(self,atom):
        # make sure cooridnates is in the deep copy
        xyz = atom.xyz[:]
        return self.Atom(atom.s,atom.n,atom.r,atom.m,atom.name,xyz)



    def pro_atom_property(self,full):
        atominfo = []
        for line in full.split('\n')[1:]:
            t = line.split(',')
            sign = t[0].strip()
            number = int(t[1])
            try:
                radius = float(t[2])
            except ValueError:
                radius = 0.0
            try:
                mass = float(t[3])
            except ValueError:
                mass = 1.0
            name = t[4].strip()
            atominfo.append([sign,number,radius,mass,name])
        return atominfo



    def get_atom(self,s=None,n=None,r=None,m=None,name=None,xyz=None,msg=True):
        if s is None:
            if n is None:
                if name is None:
                    choice = 'dummy'
                else:
                    choice = 'name'
            else:
                choice = 'number'
        else:
            choice = 'sign'

        if choice == 'dummy':
            s = 'X'
            n = 0
            r = 0.0
            m = 1.0
            name = 'Dummy'
            xyz = [0.0, 0.0, 0.0]
            return self.Atom(s,n,r,m,name,xyz)


        if choice == 'sign':
            signlist = [i[0].lower() for i in self.atominfo]
            if s.lower() in signlist:
                ndx = signlist.index(s.lower())
            else:
                numberstrlist = [str(i[1]) for i in self.atominfo]
                if s in numberstrlist:
                    ndx = numberstrlist.index(s)
                else:
                    ndx = True
                    info = 'sign < {:} >'.format(s)
        elif choice == 'number':
            numberlist = [i[1] for i in self.atominfo]
            if n in numberlist:
                ndx = numberlist.index(n)
            else:
                ndx = True
                info = 'number < {:} >'.format(n)
        else:
            namelist = [i[4].lower() for i in self.atominfo]
            if name.lower() in namelist:
                ndx = namelist.index(name.lower())
            else:
                ndx = True
                info = 'name < {:} >'.format(name)

        if ndx is True:
            if msg:
                print('Warning: not found: {:}, dummy atom is used'.format(info))
            s = 'X'
            n = 0
            r = 0.0
            m = 1.0
            name = 'Dummy'
            xyz = [0.0, 0.0, 0.0]
            return self.Atom(s,n,r,m,name,xyz)

        at = self.atominfo[ndx]
        s = at[0]
        n = at[1]
        r = at[2]
        m = at[3]
        name = at[4]

        # check xyz
        bo = False
        if xyz is None:
            xyz = [0.0, 0.0, 0.0]
        elif isinstance(xyz, list) and len(xyz) == 3:
            try:
                x = float(xyz[0])
                y = float(xyz[1])
                z = float(xyz[2])
                xyz = [x,y,z]
            except ValueError:
                bo = True
        else:
            bo = True
        if bo:
            if msg:
                print('Warning: wrong format: xyz, [0.0, 0.0, 0.0] is used')
            xyz = [0.0, 0.0, 0.0]

        return self.Atom(s,n,r,m,name,xyz)



    def guess_atoms_for_system(self,system):
        """
        Args:
            system : System[Mol[Atom, ...]] : List[List[atomtype, x,y,z], ...]]

        Return:
            sysnew : System[Mol[Atom, ...], ...]
        """
        sysnew = []
        signlist = [i[0].lower() for i in self.atominfo]
        numlist = [i[1] for i in self.atominfo]
        numstrlist = [str(i) for i in numlist]
        namelist = [i[4].lower() for i in self.atominfo]
        for mol in system:
            ls = []
            for at in mol:
                if at[0].lower in namelist:
                    ndx = namelist.index(at[0].lower())
                else:
                    if len(at[0]) <= 2:
                        tmp = at[0].lower()
                    else:
                        tmp = at[0][:2].lower()
                    if tmp in signlist:
                        ndx = signlist.index(tmp)
                    elif tmp in numlist:
                        ndx = numlist.index(tmp)
                    elif tmp in numstrlist:
                        ndx = numstrlist.index(tmp)
                    else:
                        print('Warning: not found: dummy atom will be used: {:}'.format(at[0]))
                        ndx = 0
                ls.append(self.get_atom(s=signlist[ndx],xyz=at[1:]))
            sysnew.append(ls)

        return sysnew




def test_class_AtomInfo():
    RF = ReadFileMultiple('template-mol.txt')
    RF.run()
    for mol in RF.system:
        for at in mol:
            print(at)
        print()
        print()
    ai = AtomInfo()
    system = ai.guess_atoms_for_system(RF.system)
    for mol in system:
        for at in mol:
            print(at.__dict__)
        print()
        print()




FAI = AtomInfo()




class ClosedForm:
    """Calculate best fit structure for two sets coordinates

    Args:
        vl : 2D n*3f : List[[float, float, float]] : Left input matrix
        vr : 2D n*3f : List[[float, float, float]] : Right input matrix

        vl = BestFit * vr
        Left set matrix is used for reference

    Methods:
        centroid
        calc_N
        calc_left_rotM
        calc_bestfit

    Reference:
        Closed-form solution of absolute orientation using unit quaternions
        Berthold K. P. Horn
        J. Opt. Soc. Am. A 4, 629-642 (1987)
    """
    def __init__(self,vl=None,vr=None):
        self.vl = vl
        self.vr = vr



    def calc_bestfit(self,vl=None,vr=None,centroid=True):
        """calc best fit structure

        both vl(reference) and vr(target) can be reused

        centroid: whether centroid results or not
        """
        if vl is None: vl = self.vl
        if vr is None: vr = self.vr

        cvl,tl = self.centroid(vl)
        cvr,tr = self.centroid(vr)

        M = self.calc_left_rotM(cvl,cvr)
        fit = []
        for v in cvr:
            rx = v[0]*M[0][0] + v[1]*M[1][0] + v[2]*M[2][0]
            ry = v[0]*M[0][1] + v[1]*M[1][1] + v[2]*M[2][1]
            rz = v[0]*M[0][2] + v[1]*M[1][2] + v[2]*M[2][2]
            fit.append([rx,ry,rz])

        if not centroid:
            for i in range(len(fit)):
                fit[i][0] += tr[0]
                fit[i][1] += tr[1]
                fit[i][2] += tr[2]

        return fit



    def centroid(self,v):
        """calc centroid vector and translation vector
        """
        x = [i[0] for i in v]
        y = [i[1] for i in v]
        z = [i[2] for i in v]
        t = len(v)
        ax = sum(x) / t
        ay = sum(y) / t
        az = sum(z) / t
        # care, do not use in-place operation
        cv = [[0.0 for i in range(3)] for j in range(len(v))]
        for i in range(t):
            cv[i][0] = v[i][0] - ax
            cv[i][1] = v[i][1] - ay
            cv[i][2] = v[i][2] - az
        return cv,(ax,ay,az)



    def calc_N(self,vl,vr):
        """calc 4x4 real symmetric N matrix
        """
        XxYx = 0.0
        XxYy = 0.0
        XxYz = 0.0
        XyYx = 0.0
        XyYy = 0.0
        XyYz = 0.0
        XzYx = 0.0
        XzYy = 0.0
        XzYz = 0.0
        # careful of the sequence: X-r, Y-l
        # for rotation l = Rr
        for i,p in enumerate(vl):
            XxYx += p[0] * vr[i][0]
            XxYy += p[0] * vr[i][1]
            XxYz += p[0] * vr[i][2]

            XyYx += p[1] * vr[i][0]
            XyYy += p[1] * vr[i][1]
            XyYz += p[1] * vr[i][2]

            XzYx += p[2] * vr[i][0]
            XzYy += p[2] * vr[i][1]
            XzYz += p[2] * vr[i][2]

        N = [[0.0, 0.0, 0.0, 0.0] for i in range(4)]

        N[0][0] = XxYx + XyYy + XzYz
        N[0][1] = XyYz - XzYy
        N[0][2] = XzYx - XxYz
        N[0][3] = XxYy - XyYx

        N[1][0] = N[0][1]
        N[1][1] = XxYx - XyYy - XzYz
        N[1][2] = XxYy + XyYx
        N[1][3] = XzYx + XxYz

        N[2][0] = N[0][2]
        N[2][1] = N[1][2]
        N[2][2] = -XxYx + XyYy - XzYz
        N[2][3] = XyYz + XzYy

        N[3][0] = N[0][3]
        N[3][1] = N[1][3]
        N[3][2] = N[2][3]
        N[3][3] = -XxYx - XyYy + XzYz

        return N



    def calc_left_rotM(self,vl,vr):
        """calc rotation matrix for vr*M computed from quaternion

        quaternion is got from the vector corresponding to
        largest positive eigenvalue

        M  : 2D 3*3f : rotation matrix for vr, vr*M, in element-wise operation
        """
        N = self.calc_N(vl,vr)
        values,vectors = np.linalg.eig(N)
        ndx = np.where(values == max(values))
        ndx = ndx[0][0]
        # For numpy, eigenvectors are correspondingly put in column
        # note, this vector has already been normalized
        V = vectors[:,ndx]
        M = [[0.0, 0.0, 0.0] for i in range(3)]

        M[0][0] = 1 - 2 * (V[2]*V[2] + V[3]*V[3])
        M[0][1] = 2 * (V[1]*V[2] - V[3]*V[0])
        M[0][2] = 2 * (V[0]*V[3] + V[2]*V[0])

        M[1][0] = 2 * (V[1]*V[2] + V[3]*V[0])
        M[1][1] = 1 - 2 * (V[1]*V[1] + V[3]*V[3])
        M[1][2] = 2 * (V[2]*V[3] - V[1]*V[0])

        M[2][0] = 2 * (V[1]*V[3] - V[2]*V[0])
        M[2][1] = 2 * (V[2]*V[3] + V[1]*V[0])
        M[2][2] = 1 - 2 * (V[1]*V[1] + V[2]*V[2])

        return M




def test_class_ClosedForm():
    CF = ClosedForm()

    vl = [[0,0,0],[0,1,0],[1,0,0]]
    vr = [[0,0,0],[0,1,0],[0,0,1]]
    rst = CF.calc_bestfit(vl,vr)
    cvl,t = CF.centroid(vl)
    assert np.allclose(cvl,rst)

    vl = [[0,0,0],[0,1,0]]
    vr = [[0,0,0],[0,0,1]]
    rst = CF.calc_bestfit(vl,vr)
    cvl,t = CF.centroid(vl)
    assert np.allclose(cvl,rst)

    vl = [[0,0,0],[0,1,1]]
    vr = [[0,0,0],[0,0,pow(2,0.5)]]
    rst = CF.calc_bestfit(vl,vr)
    cvl,t = CF.centroid(vl)
    assert np.allclose(cvl,rst)

    vl = [[0,0,0],[0,1,0],[1,0,0],[0,1,0]]
    vr = [[0,0,0],[0,1,0],[0,0,1],[0,1,0]]
    rst = CF.calc_bestfit(vl,vr)
    cvl,t = CF.centroid(vl)
    assert np.allclose(cvl,rst)

    vl = [[0,0,0],[0,1,0],[1,0,0],[0,1,0],[1,1,0]]
    vr = [[0,0,0],[0,1,0],[0,0,1],[0,1,0],[0,1,1]]
    rst = CF.calc_bestfit(vl,vr)
    cvl,t = CF.centroid(vl)
    assert np.allclose(cvl,rst)




def rotation(pd=None,po=None,angle=None,v=None):
    """rotation matrix for axis po->pd in given angle

    Args:
        pd : 1D 1*3f, direction point, default [0,0,1]
        po : 1D 1*3f, origin point, default [0,0,0]
        angle (float): in degree, default 5.0

        v  : 1D 1*3f, vector, it should be proportional to po->pd, low priority

    Return:
        R  : 1D 4*4f : built by [Rreal, T], thus targets format [[x,y,z], 1], matrix operation

    Reference:
        Rotation About an Arbitrary Axis in 3 Dimensions
        Glenn Murray
        https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas
        /rotation-about-an-arbitrary-axis-in-3-dimensions
    """
    if po is None and pd is None and v is not None:
        po = [0.0, 0.0, 0.0]
        pd = v
    else:
        po = [0.0, 0.0, 0.0] if po is None else po
        pd = [0.0, 0.0, 1.0] if pd is None else pd
    angle = 5.0 if angle is None else angle

    angle = angle / 180.0 * np.pi
    l = [pd[i]-po[i] for i in range(3)]
    tt = sum([i*i for i in l])
    t = pow(tt,0.5)
    u,v,w = l[0]/t, l[1]/t, l[2]/t
    a,b,c = po[0],po[1],po[2]
    cg = np.cos(angle)
    sg = np.sin(angle)

    R = [[0, 0, 0, 0] for i in range(4)]
    R[0][0] = u*u + (v*v + w*w)*cg
    R[0][1] = u*v*(1-cg) - w*sg
    R[0][2] = u*w*(1-cg) + v*sg
    R[0][3] = (a*(v*v+w*w) - u*(b*v+c*w))*(1-cg) + (b*w-c*v)*sg

    R[1][0] = u*v*(1-cg) + w*sg
    R[1][1] = v*v + (u*u+w*w)*cg
    R[1][2] = v*w*(1-cg) - u*sg
    R[1][3] = (b*(u*u+w*w) - v*(a*u+c*w))*(1-cg) + (c*u-a*w)*sg

    R[2][0] = u*w*(1-cg) - v*sg
    R[2][1] = v*w*(1-cg) + u*sg
    R[2][2] = w*w + (u*u+v*v)*cg
    R[2][3] = (c*(u*u+v*v) - w*(a*u+b*v))*(1-cg) + (a*v-b*u)*sg

    R[3][0] = 0.0
    R[3][1] = 0.0
    R[3][2] = 0.0
    R[3][3] = 1.0

    return R




def test_func_rotation():
    def multi(R,p):
        x = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2] + R[0][3]
        y = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2] + R[1][3]
        z = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2] + R[2][3]
        return [x,y,z]

    pd = [0,0,1]
    pt = [1,0,0]
    angle = 90
    expect = [0,1,0]
    R = rotation(pd=pd,angle=angle)
    assert np.allclose(expect,multi(R,pt))

    pd = [0,0,1]
    pt = [pow(2,0.5),0,0]
    angle = 135
    expect = [-1,1,0]
    R = rotation(pd=pd,angle=angle)
    assert np.allclose(expect,multi(R,pt))

    pd = [1,1,0]
    pt = [1,-1,0]
    angle = 90
    expect = [0,0,-pow(2,0.5)]
    R = rotation(pd=pd,angle=angle)
    assert np.allclose(expect,multi(R,pt))

    pd = [0,1,0]
    pt = [0,0,1]
    angle = 90
    expect = [1,0,0]
    R = rotation(pd=pd,angle=angle)
    assert np.allclose(expect,multi(R,pt))

    po = [0,0,0.5]
    pd = [0,1,0.5]
    pt = [0,0,1]
    angle = 90
    expect = [0.5,0,0.5]
    R = rotation(pd=pd,po=po,angle=angle)
    assert np.allclose(expect,multi(R,pt))

    po = [0,0,-1]
    pd = [1,0,-1]
    pt = [0,0,1]
    angle = 60
    expect = [0,-pow(3,0.5),0]
    ER = [[1,0,0,0],[0,0.5,-0.866,-0.866],[0,0.866,0.5,-0.5],[0,0,0,1]]
    R = rotation(pd=pd,po=po,angle=angle)
    assert np.allclose(expect,multi(R,pt))
    assert np.allclose(R,ER,rtol=0.0001,atol=0.001)

    po = [1,7,3]
    pd = [3,16,9]
    pt = [5,8,4]
    angle = 26
    expect = [4.7532,8.9487,2.6593]
    ER = [
        [0.9021, -0.2241,  0.3687,  0.5601],
        [0.2542,  0.9665, -0.0345,  0.0836],
        [-0.3486, 0.1249,  0.9289, -0.3122],
        [0, 0, 0, 1]
    ]
    R = rotation(pd=pd,po=po,angle=angle)
    assert np.allclose(expect,multi(R,pt),rtol=0.00001,atol=0.0001)
    assert np.allclose(R,ER,rtol=0.0001,atol=0.0001)




def translation(pd=None,po=None,v=None):
    """translation matrix for vector determined by po->pd

    v  : 1D 1*3f, vector, it should be proportional to po->pd, low priority

    Return:
        T : 1D 1*3f     :   element-wise, [x, y, z]
    """
    if pd is None and po is None and v is not None:
        po = [0.0, 0.0, 0.0]
        pd = v
    else:
        po = [0.0, 0.0, 0.0] if po is None else po
        pd = [1.0, 0.0, 0.0] if pd is None else pd

    tol = 0.0000000001

    # calculate their sin(angle) values
    v = [pd[i] - po[i] for i in range(3)]
    vv = [i*i for i in v]

    xyz = sum(vv)
    xy = sum(vv[:2])

    sxy = pow(vv[2]/xyz, 0.5) if xyz > tol else tol
    syz = pow(vv[0]/xy, 0.5) if xy > tol else tol
    sxz = pow(vv[1]/xy, 0.5) if xy > tol else tol

    if v[0] < 0: syz = -syz
    if v[1] < 0: sxz = -sxz
    if v[2] < 0: sxy = -sxy

    return [syz, sxz, sxy]




def test_func_translation():
    pt = [1,1,0]

    T = translation()
    ptsnew = []
    for i in range(10):
        nx = pt[0] + T[0]*i
        ny = pt[1] + T[1]*i
        nz = pt[2] + T[2]*i
        ptsnew.append([nx,ny,nz])
    for i in ptsnew:
        print(i)

    print('\n\n')
    pd = [-1,0,0]
    T = translation(pd=pd)
    ptsnew = []
    for i in range(10):
        nx = pt[0] + T[0]*i
        ny = pt[1] + T[1]*i
        nz = pt[2] + T[2]*i
        ptsnew.append([nx,ny,nz])
    for i in ptsnew:
        print(i)

    print('\n\n')
    inc = pow(2,0.5)
    pd = [-1,1,0]
    T = translation(pd=pd)
    ptsnew = []
    for i in range(10):
        nx = pt[0] + T[0]*i*inc
        ny = pt[1] + T[1]*i*inc
        nz = pt[2] + T[2]*i*inc
        ptsnew.append([nx,ny,nz])
    for i in ptsnew:
        print(i)




class ParConfig:
    """check and set parameters

    Args:
        system  : 2D n*nAtom : [ [Atom, ...],  ...]

        po   |  pd   |  pr   |  pt      :   1D 1*3f
        ipo  |  ipd  |  ipr  |  ipt     :   int

        pts     :   2D 2*3f
        ipts    :   1D 1*ni  :  nm-atoms, start at 0

        idpo : 1D 1*2i : [mol-index, atom-index], start at 0
        idpd : 1D 1*2i : [mol-index, atom-index]
        idpr : 1D 1*2i : [mol-index, atom-index]
        idpt : 1D 1*2i : [mol-index, atom-index]
        idpts : 2D n*2i: [[mol-index, atom-index], ]  : atoms to be sampled


    assume right-hand cooridnates system

         z   pt(VT)
         |  /
         | /       angle for plane: P(po-pd-pr) and P(po-pd-pt) <= 180.0
       po|/----y   pt is on the right of pr, if vector po->pd is in front/top
         /\ 
        /   \ 
    x(pd)  pr(VR)

    Attributes:
        type (str)  : bond | angle | dihedral | zoom | None

            bond    :   po -> pd
                inc = 0.01  (Angstrom)
                ratio = 2 | 0.1

            angle   :   pd -> po -> pt
                mode :  horizontal | vertical
                inc = 0.01  (degree)
                ratio = 1.8 | 0.1

            dihedral:   pr -> po -> pd -> pt
                inc = 0.01  (degree)
                ratio = 1.5 | 0.1

            zoom    :   po
                inc = 0.001  (Angstrom)
                ratio = 1.5 | 0.3

        v : 1D 1*3f : rotation axis, self.po->self.pd
        value (float): distance | angle
        num (int) : estimated number of generations

    Note:
        1) please use [type,want,mode] to configure parameters
            i) bond:
                more:
                    suggest to use ratio (>1.0)
                less:
                    suggest to use ratio (<1.0)
                    set flip = True/False, whether allows to pass its origin
            ii) angle:
                more:
                    suggest to use ratio (>1.0)
                less:
                    horizontal:
                        suggest to use ratio (>1.0)
                        con shape, no influence at all
                    vertical:
                        suggest to use ratio (<1.0)
                        set flip = True/False, whether allows to pass its origin
            iii) dihedral:
                more:
                    suggest to use ratio (>1.0)
                less:
                    suggest to use ratio (<1.0)
                    set flip = True/False, whether allows to pass its origin
            iiii) zoom:
                more:
                    suggest to use ratio (>1.0)
                less:
                    suggest to use ratio (<1.0)
                    flip is not allowed

        2) please check self.nice before using.
            *) True means no errors
            *) Otherwise, check self.info

        3) priority:
            => self.v 
                => self.po -> self.pd
                    => self.ipo -> self.ipd
                        => self.idpo -> self.idpd

        4) self.num for double-check
    """
    def __init__(self,system,*args,v=None,
                po=None,pd=None,pr=None,pt=None,pts=None,
                ipo=None,ipd=None,ipr=None,ipt=None,ipts=None,
                idpo=None,idpd=None,idpr=None,idpt=None,idpts=None,
                type=None,want=None,mode=None,flip=None,
                inc=None,ratio=None,threshold=None,**kwargs):

        self.nice = True
        self.info = ''

        self.system = system
        if self.system is None:
            self.nice = False
            self.info = 'Fatal: not defined: system'
            return
        # calculate its atoms accumulation list
        tmplist = [len(i) for i in self.system]
        self._acclist = [0,]
        tot = 0
        for i in tmplist:
            tot += i
            self._acclist.append(tot)

        self.v = v
        if self.v is not None:
            if not isinstance(v,list) or len(v) != 3:
                self.nice = False
                self.info = 'Fatal: wrong defined: v: {:}'.format(v)
                return

        if self._check('po',po): return
        if self._check('pd',pd): return
        if self._check('pr',pr): return
        if self._check('pt',pt): return
        if pts is not None:
            for p in pts:
                if self._check('pts',p): return

        if self._check_i('ipo',ipo): return
        if self._check_i('ipd',ipd): return
        if self._check_i('ipr',ipr): return
        if self._check_i('ipt',ipt): return
        if ipts is not None:
            for p in ipts:
                if self._check_i('ipts',p): return

        if self._check_id('idpo',idpo): return
        if self._check_id('idpd',idpd): return
        if self._check_id('idpr',idpr): return
        if self._check_id('idpt',idpt): return
        if idpts is not None:
            for p in idpts:
                if self._check_id('idpts',p): return

        self.po = po
        self.pd = pd
        self.pr = pr
        self.pt = pt
        self.pts = pts if pts is not None else []
        
        self.ipo = ipo
        self.ipd = ipd
        self.ipr = ipr
        self.ipt = ipt
        self.ipts = ipts

        self.idpo = idpo
        self.idpd = idpd
        self.idpr = idpr
        self.idpt = idpt
        self.idpts = idpts

        # conclude all target points into idpts
        if idpts is None: self.idpts = []
        if self.ipts is not None: self.idpts.extend(self.calc_pid(self.ipts))
        if len(self.idpts) + len(self.pts) == 0:
            self.nice = False
            self.info = 'Fatal: not defined: not target points'
            return
        if len(self.idpts) != 0:
            # remove repeats
            reflist = []
            for i,p in enumerate(self.idpts):
                bo = True
                for j in range(i+1,len(self.idpts)):
                    if self.idpts[j][0] == p[0] and self.idpts[j][1] == p[1]:
                        bo = False
                        break
                if bo:
                    reflist.append(i)
            self.idpts = [self.idpts[i] for i in reflist]

        self.type = type
        if self.type is None:
            if pr is None:
                if pt is None:
                    self.type = 'bond'
                else:
                    self.type = 'angle'
            else:
                self.type = 'dihedral'

        self.mode = mode
        self.flip = True if flip is True else False
        self.inc = 0.01 if inc is None else max(inc, 0.000001) # avoid zero division
        
        if ratio is not None and not isinstance(ratio,(float,int)):
            self.nice = False
            self.info = 'Fatal: not a number: ratio: {:}'.format(ratio)
            return
        self.ratio = ratio

        if threshold is not None and not isinstance(threshold,(float,int)):
            self.nice = False
            self.info = 'Fatal: not a number: threshold: {:}'.format(threshold)
            return
        self.threshold = threshold

        if want is None:
            if self.ratio is None:
                self.want = 'more'
            elif self.ratio > 1.0:
                self.want = 'more'
            else:
                self.want = 'less'
        elif want.lower() in ['more','bigger','longer','larger','anticlockwise','m','b','a']:
            self.want = 'more'
        elif want.lower() in ['less','smaller','shorter','clockwise','l','s','c']:
            self.want = 'less'
        else:
            self.nice = False
            self.info = 'Fatal: wrong configure: want: {:}'.format(want)
            return

        # set po
        if self.po is None:
            if self.ipo is None:
                if self.idpo is None:
                    self.po = [0.0, 0.0, 0.0]
                else:
                    self.po = self.system[self.idpo[0]][self.idpo[1]].xyz
            else:
                self.idpo = self.calc_pid(ipo)
                self.po = self.system[self.idpo[0]][self.idpo[1]].xyz

        # set pd
        if self.pd is None:
            if self.ipd is None:
                if self.idpd is None:
                    self.pd = [1.0, 0.0, 0.0]
                else:
                    self.pd = self.system[self.idpd[0]][self.idpd[1]].xyz
            else:
                self.idpd = self.calc_pid(ipd)
                self.pd = self.system[self.idpd[0]][self.idpd[1]].xyz

        # set pr
        if self.pr is None:
            if self.ipr is None:
                if self.idpr is None:
                    self.pr = [0.0, 1.0, 0.0]
                else:
                    self.pr = self.system[self.idpr[0]][self.idpr[1]].xyz
            else:
                self.idpr = self.calc_pid(ipr)
                self.pr = self.system[self.idpr[0]][self.idpr[1]].xyz

        # set pt
        if self.pt is None:
            if self.ipt is None:
                if self.idpt is None:
                    self.pt = [0.0, 0.0, 1.0]
                else:
                    self.pt = self.system[self.idpt[0]][self.idpt[1]].xyz
            else:
                self.idpt = self.calc_pid(ipt)
                self.pt = self.system[self.idpt[0]][self.idpt[1]].xyz

        # set v
        if self.v is None:
            self.v = [self.pd[i]-self.po[i] for i in range(3)]
        else:
            self.po = [0.0, 0.0, 0.0]
            self.pd = [i for i in self.v]

        if self.type.lower() in ['bond','bonds','b']:
            self.type = 'bond'
        elif self.type.lower() in ['angle','angles','a']:
            self.type = 'angle'
        elif self.type.lower() in ['dihedral','dihedrals','d']:
            self.type = 'dihedral'
        elif self.type.lower() in ['zoom','zooms','z']:
            self.type = 'zoom'
            # care change zoom default inc
            if inc is None: self.inc = 0.001
        else:
            self.info = 'Warning: wrong configure: type: {:}'.format(self.type)



    def _check(self,pid,p):
        if p is None:
            return False
        if not isinstance(p,list):
            self.nice = False
            self.info = 'Fatal: not a list: {:}: {:}'.format(pid,p)
            return True
        if len(p) != 3:
            self.nice = False
            self.info = 'Fatal: length should be 3 [x,y,z]: {:}: {:}'.format(pid,p)
            return True
        for i in p:
            if not isinstance(i,(float,int)):
                self.nice = False
                self.info = 'Fatal: not a number: {:}: {:}'.format(pid,p)
                return True
        return False



    def _check_i(self,pid,p):
        if p is None:
            return False
        if not isinstance(p,int) or p < 0:
            self.nice = False
            self.info = 'Fatal: not an integer: {:}: {:}'.format(pid,p)
            return True
        return False



    def _check_id(self,pid,p):
        if p is None:
            return False
        if not isinstance(p,list):
            self.nice = False
            self.info = 'Fatal: not a list: {:}: {:}'.format(pid,p)
            return True
        if len(p) != 2:
            self.nice = False
            self.info = 'Fatal: length should be 2 [mol-ndx, atom-ndx]: {:}: {:}'.format(pid,p)
            return True
        if len(self.system) < p[0]:
            self.nice = False
            self.info = 'Fatal: molecule index too large: {:}: {:}'.format(pid,p)
            return True
        if len(self.system[p[0]]) < p[1]:
            self.nice = False
            self.info = 'Fatal: atom index too large: {:}: {:}'.format(pid,p)
            return True
        return False



    def config_bond(self):
        dd = sum([i*i for i in self.v])
        self.value = pow(dd,0.5)

        if self.want == 'more':
            self.lower = self.value
            self.info = 'bond longer/clockwise'
            if self.threshold is None:
                tmp = 2.0 if self.ratio is None else self.ratio
                self.higher = self.value * tmp
            else:
                self.higher = self.threshold
        else:
            self.higher = self.value
            self.info = 'bond shorter/anticlockwise'
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]
            if self.threshold is None:
                tmp = 0.1 if self.ratio is None else self.ratio
                if self.flip is True: tmp = -tmp
                self.lower = self.value * tmp
            else:
                if self.flip is True: self.threshold = -self.threshold
                self.lower = self.threshold

        self._calc_num()



    def config_angle(self):
        """
        Given three points, po,pd,pt, to make angle pd->po->pt right-handed,
        which means, assume angle less than 180 degree, direction point
        pd is placed on the left of target point pt, if they are in oxy planar

        mathematically:
            pt(y) > l(x)    where l is the function for line po->pd


        suppose two vectors (v1, v2, v3), (s1, s2, s3),
        if we want to find another vector (x, y, z) that is perpendicular to
        those two vectors, we will have:

               z   pt(S)
            V  |  /
               | /
               o ----- y      points o,pt,pd are inside planar oxy
              / \             this is so-called righthand coordinates system
             /   \ 
            x     pd(H)

                [  x   y   z   ]            x             y             z
        H x S = |  h1  h2  h3  |  =  ( h2*s3-h3*s2   h3*s1-h1*s3   h1*s2-h2*s1 )
                [  s1  s2  s3  ]
        """
        yc = self.v[1] * (self.pt[0]-self.po[0])
        yr = self.v[0] * (self.pt[1]-self.po[1])
        if yc > yr:
            # update vector po->pd to pd->po
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]


        if self.mode in ['vertical','v']:
            bo = True
        elif self.mode in ['horizontal','h']:
            bo = False
        else:
            # redundant codes for future updates
            bo = False

        s = [self.pt[i]-self.po[i] for i in range(3)]
        tt = sum([i*i for i in s])
        dd = sum([i*i for i in self.v])
        dt = sum([self.v[i]*s[i] for i in range(3)])
        tmp = max(pow(dd*tt,0.5), 0.000000000001)     # to avoid zero division
        self.value = np.arccos(dt/tmp) * 180.0 / np.pi
        if bo:
            self.mode = 'vertical'
            x = self.v[1]*s[2] - self.v[2]*s[1]
            y = self.v[2]*s[0] - self.v[0]*s[2]
            z = self.v[0]*s[1] - self.v[1]*s[0]
            self.v = [x, y, z]
        else:
            self.mode = 'horizontal'

        self.lower = self.value
        if self.threshold is None:
            tmp = 1.8 if self.ratio is None else self.ratio
            self.higher = self.value * tmp
        else:
            self.higher = self.threshold
        if self.higher < self.lower:
                self.lower, self.higher = self.higher, self.lower
        
        self.info = 'angle rotation, righthandness, more/anticlockwise'
        if self.mode == 'horizontal':
            if self.want == 'less':
                self.info = 'angle rotation, righthandness, less/clockwise'
                self.po, self.pd = self.pd, self.po
                self.v = [-i for i in self.v]
        else:
            if self.want == 'less':
                self.info = 'angle rotation, righthandness, less/clockwise'
                self.higher = self.value
                if self.threshold is None:
                    tmp = 0.1 if self.ratio is None else self.ratio
                    self.lower = self.value * tmp
                else:
                    self.lower = self.threshold
                if self.flip: self.lower = -self.lower
                self.po, self.pd = self.pd, self.po
                self.v = [-i for i in self.v]

        self._calc_num()



    def config_dihedral(self):
        """
        Similarly for angle configuration, set points pr -> vector(po-pd) -> pt
        are in right-handed system. For example,

            z  pt(VT)
            | /
            |/ po     angle for plane: P(po-pd-pr) and P(po-pd-pt) <= 180.0
            |/----y   pt is on the right of pr, if vector po->pd faces to front/top
           /\ 
          /   \ 
        x(pd)  pr(VR)

        hard to explain, value can be larger than 180.0
        """
        # first, set pd->po->pr right-handed
        yc = self.v[1] * (self.pr[0]-self.po[0])
        yr = self.v[0] * (self.pr[1]-self.po[1])
        if yc > yr:
            # update vector po->pd to pd->po
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]

        # second, calculate norm vector for plane pd->po->pr
        vr = [self.pr[i]-self.po[i] for i in range(3)]
        x = self.v[1]*vr[2] - self.v[2]*vr[1]
        y = self.v[2]*vr[0] - self.v[0]*vr[2]
        z = self.v[0]*vr[1] - self.v[1]*vr[0]
        nr = [x, y, z]
        rr = sum([i*i for i in nr])

        # third, calculate norm vector for plane pd->po->pt, be aware!
        # only in this way, the dihedral will be equal to vector angle
        vt = [self.pt[i]-self.po[i] for i in range(3)]
        x = self.v[1]*vt[2] - self.v[2]*vt[1]
        y = self.v[2]*vt[0] - self.v[0]*vt[2]
        z = self.v[0]*vt[1] - self.v[1]*vt[0]
        nt = [x, y, z]
        tt = sum([i*i for i in nt])

        # fourth, calculate dihedral
        tr = sum([nt[i]*nr[i] for i in range(3)])
        tmp = max(pow(tt*rr,0.5), 0.000000000001)     # to avoid zero division
        self.value = np.arccos(tr/tmp) * 180.0 / np.pi

        # fifth, make pt->v->pr right-handed, where v is vector po->pd or pd->po,
        # which means, for vector v, pt is right-handed for pr, considering
        # dihedral between plane po,pd,pr and plane po,pd,pt is less than 180.0
        zc = nr[2] * nt[1]
        zr = nr[1] * nt[2]
        if zc > zr:
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]

        if self.want == 'more':
            self.info = 'dihedral rotation, righthandness, more/anticlockwise'
            self.lower = self.value
            if self.threshold is None:
                tmp = 1.5 if self.ratio is None else self.ratio
                self.higher = self.value * tmp
            else:
                self.higher = self.threshold
        else:
            self.info = 'dihedral rotation, righthandness, less/clockwise'
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]
            self.higher = self.value
            if self.threshold is None:
                tmp = 0.1 if self.ratio is None else self.ratio
                if self.flip is True: tmp = -tmp
                self.lower = self.value * tmp
            else:
                if self.flip is True: self.threshold = -self.threshold
                self.lower = self.threshold

        self._calc_num()



    def config_zoom(self):
        dd = sum([i*i for i in self.v])
        self.value = pow(dd,0.5)
        self.flip = False

        if self.want == 'more':
            self.lower = self.value
            self.info = 'zoom in'
            if self.threshold is None:
                tmp = 1.5 if self.ratio is None else self.ratio
                self.higher = self.value * tmp
            else:
                self.higher = self.threshold
        else:
            self.higher = self.value
            self.info = 'zoom out'
            self.po, self.pd = self.pd, self.po
            self.v = [-i for i in self.v]
            if self.threshold is None:
                tmp = 0.3 if self.ratio is None else self.ratio
                self.lower = self.value * tmp
            else:
                self.lower = self.threshold

        self._calc_num()



    def _calc_num(self):
        if self.lower >= self.higher:
            self.num = 0
            self.info += '\nWarning: wrong setting: threshold/ratio: no generation'
            return
        self.num = (self.higher-self.lower) / self.inc
        self.num = int(self.num+0.51)
    


    def calc_pid(self,pis):
        """calculate integer to index to [mol-index, atom-index]

        Args:
            pis: int | 1D 1*ni  :   counting starts at 0

        Return:
            1D 1*2i     [int, int]          for int
            2D n*2i     [[int,int], ...]    Otherwise
        """
        if isinstance(pis,int):
            for i,nm in enumerate(self._acclist):
                if pis < nm:
                    mi = i - 1
                    ai = pis - self._acclist[i-1]
                    break
            return [mi,ai]

        reflist = []
        for atnum in pis:
            for i,nm in enumerate(self._acclist):
                if atnum < nm:
                    mi = i - 1
                    ai = atnum - self._acclist[i-1]
                    break
            reflist.append([mi,ai])
        return reflist




def test_class_ParConfig():
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000

        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    rf = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    rf.run()

    system = FAI.guess_atoms_for_system(rf.system)
    for mol in system:
        for at in mol:
            print(at.__dict__)
        print('\n')


    # test bond default
    pts = [[9,0,0]]
    pc = ParConfig(system,pts=pts)
    assert pc.nice
    assert pc.type == 'bond'
    pc.config_bond()
    assert pc.want == 'more'

    # test bond i -> id
    idpo = [1,1]
    idpd = [0,1]
    idpts = [[0,1],[1,1]]
    pc = ParConfig(system,idpo=idpo,idpd=idpd,pts=pts,ipts=[1,4])
    assert pc.nice
    assert pc.type == 'bond'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert np.allclose(pc.po, system[1][1].xyz)
    assert np.allclose(pc.pd, system[0][1].xyz)
    assert np.allclose(pc.idpts, idpts)

    # test bond ratio, switch po <-> pd
    pc = ParConfig(system,ratio=0.3,idpo=idpo,ipd=2,idpd=idpd,pts=pts,ipts=[1,4])
    assert pc.nice
    assert pc.type == 'bond'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'less'
    assert np.allclose(pc.po, system[0][2].xyz)
    assert np.allclose(pc.pd, system[1][1].xyz)

    # test bond priority, v > pd > ipd > idpd
    pc = ParConfig(system,v=[1,0,0],pd=[8,0,0],idpo=idpo,ipd=2,idpd=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'bond'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert np.allclose(pc.v, [1,0,0])

    pc = ParConfig(system,pd=[8,0,0],idpo=idpo,ipd=2,idpd=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'bond'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert np.allclose(pc.pd, [8,0,0])


    # test angle default
    pc = ParConfig(system,pt=[1,1,1],pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.mode == 'horizontal'
    assert np.allclose(pc.v, [1,0,0])

    pc = ParConfig(system,pt=[0,1,0],ratio=0.3,mode='v',flip=True,pts=pts)
    assert pc.nice
    getattr(pc,'config_'+pc.type)()
    assert pc.flip
    assert pc.mode == 'vertical'
    assert pc.want == 'less'
    assert np.allclose(pc.v, [0,0,-1])

    # test angle priority, v > pt > ipt > idpt
    pc = ParConfig(system,type='a',v=[7,0,0],idpo=idpo,ipt=2,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert np.allclose(pc.v, [7,0,0])

    pc = ParConfig(system,type='a',v=[-7,0,0],flip=True,idpo=idpo,ipt=1,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.flip
    assert np.allclose(pc.pt, system[0][1].xyz)
    assert np.allclose(pc.v, [7,0,0])

    pc = ParConfig(system,type='a',mode='v',pt=[8,0,0],idpo=idpo,ipt=2,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.mode == 'vertical'
    assert np.allclose(pc.pt, [8,0,0])

    pc = ParConfig(system,type='a',idpo=idpo,ipt=2,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.mode == 'horizontal'
    assert np.allclose(pc.pt, system[0][2].xyz)

    pc = ParConfig(system,type='a',mode='v',v=[1,1,1],idpo=idpo,ipt=1,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.mode == 'vertical'
    assert pc.v[2] > 0

    pc = ParConfig(system,type='a',mode='v',v=[-1,-1,1],idpo=idpo,ipt=1,idpt=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'angle'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert pc.mode == 'vertical'
    assert pc.v[2] > 0
 

    # test dihedral default
    pc = ParConfig(system,pr=[6,0,0],pts=pts)
    assert pc.nice
    assert pc.type == 'dihedral'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'

    pc = ParConfig(system,ratio=0.1,pt=[8,0,0],pr=[7,0,0],pts=pts)
    assert pc.nice
    assert pc.type == 'dihedral'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'less'

    # test dihedral priority, v > pr > ipr > idpr
    pc = ParConfig(system,type='d',v=[7,0,0],idpo=idpo,ipr=1,idpr=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'dihedral'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'more'
    assert np.allclose(pc.v, [7,0,0])
    assert np.allclose(pc.pr, system[0][1].xyz)

    pc = ParConfig(system,type='d',want='less',v=[7,0,0],idpo=idpo,ipr=1,idpr=idpd,pts=pts)
    assert pc.nice
    assert pc.type == 'dihedral'
    getattr(pc,'config_'+pc.type)()
    assert pc.want == 'less'
    assert np.allclose(pc.v, [-7,0,0])
    assert np.allclose(pc.pr, system[0][1].xyz)




class VaryBond(ParConfig):
    """Bonds variation
    """
    def __init__(self,system,*args,**kwargs):
        super().__init__(system,*args,**kwargs)
        if not self.nice: return
        if self.type != 'bond':
            self.nice = False
            self.info += '\nFatal: configure error: not bonds variation'
            return
        getattr(self,'config_'+self.type)()



    def run(self):
        self.sysnew = []
        self.syspts = []
        mat = translation(v=self.v)
        for i in range(1,self.num+1):
            if len(self.idpts) != 0:
                mol = [[FAI.copy_atom(i) for i in j] for j in self.system]
                for ndx in self.idpts:
                    mol[ndx[0]][ndx[1]].xyz[0] += self.inc * i * mat[0]
                    mol[ndx[0]][ndx[1]].xyz[1] += self.inc * i * mat[1]
                    mol[ndx[0]][ndx[1]].xyz[2] += self.inc * i * mat[2]
                self.sysnew.append(mol)

            if len(self.pts) != 0:
                ls = []
                for pt in self.pts:
                    x = pt[0] + self.inc * i * mat[0]
                    y = pt[1] + self.inc * i * mat[1]
                    z = pt[2] + self.inc * i * mat[2]
                    ls.append([x,y,z])
                self.syspts.append(ls)





def test_class_VaryBond():
    system = [[FAI.get_atom(s='h'), FAI.get_atom(s='c',xyz=[1,1,1])],
            [FAI.get_atom(s='o',xyz=[1,2,3]), FAI.get_atom(s='s',xyz=[0,1,1])]]

    for i in system:
        for j in i:
            print(j.__dict__)
        print()

    idpts = [[0,1],[1,1]]

    fv = VaryBond(system,idpts=idpts,pd=[0,0,1])
    if not fv.nice:
        print(fv.info)
        return
    fv.run()
    print(fv.idpts)

    print(fv.num)
    for s in range(min(fv.num,5)):
        for m in fv.sysnew[s]:
            for a in m:
                print(a.__dict__)
        print()
        print()




class VaryAngle(ParConfig):
    """Angles variation
    """
    def __init__(self,system,*args,**kwargs):
        super().__init__(system,*args,**kwargs)
        if not self.nice: return
        if self.type != 'angle':
            self.nice = False
            self.info += '\nFatal: configure error: not angles variation'
            return
        getattr(self,'config_'+self.type)()



    def run(self):
        self.sysnew = []
        self.syspts = []
        for i in range(1,self.num+1):
            if len(self.idpts) != 0:
                mol = [[FAI.copy_atom(i) for i in j] for j in self.system]
                mat = rotation(v=self.v,angle=self.inc*i)
                for ndx in self.idpts:
                    mol[ndx[0]][ndx[1]].xyz = self.quaternion(mol[ndx[0]][ndx[1]].xyz, mat)
                self.sysnew.append(mol)

            if len(self.pts) != 0:
                ls = []
                for pt in self.pts:
                    ls.append( self.quaternion(pt, mat) )
                self.syspts.append(ls)



    def quaternion(self,p,R):
        """
        Args:
            p   : 1D 1*3f
            R   : 2D 4*4f   : rotation matrix

        Return:
            vls : 2D n*3f
        """
        x = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2] + R[0][3]
        y = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2] + R[1][3]
        z = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2] + R[2][3]
        return [x, y, z]



def test_class_VaryAngle():
    system = [[FAI.get_atom(s='h'), FAI.get_atom(s='c',xyz=[1,1,3])],
            [FAI.get_atom(s='o',xyz=[0,1,1]), FAI.get_atom(s='s',xyz=[0,1,1])]]

    pt = [1,1,1]
    for i in system:
        for j in i:
            print(j.__dict__)
        print()

    idpts = [[0,1],[1,1]]

    fv = VaryAngle(system,type='a',idpts=idpts,ratio=3.0)
    print(fv.v)
    print(fv.want)
    print(fv.idpts)
    if not fv.nice:
        print(fv.info)
        return
    fv.run()
    print(fv.lower,fv.higher)
    print(fv.po, fv.pt, fv.pd)

    ls = []
    for i in range(fv.num):
        if i % 300 == 0:
            ls.append(fv.sysnew[i][0][1].xyz)

    for i in range(min(len(ls),20)):
        print(ls[i],'===>   ', sum([j*j for j in ls[i]]))

    print(fv.num)
    xl = []
    yl = []
    zl = []
    for i in range(fv.num):
        if i % 300 == 0:
            for mol in fv.sysnew[i]:
                for at in mol:
                    xl.append(at.xyz[0])
                    yl.append(at.xyz[1])
                    zl.append(at.xyz[2])

    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(xl,yl,zl)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    plt.close()




class VaryDihedral(VaryAngle):
    """Dihedrals variation
    """
    def __init__(self,system,*args,**kwargs):
        # be aware in here, we are calling its grandparent
        ParConfig.__init__(self,system,*args,**kwargs)
        if not self.nice: return
        if self.type != 'dihedral':
            self.nice = False
            self.info += '\nFatal: configure error: not dihedrals variation'
            return
        getattr(self,'config_'+self.type)()




def test_class_VaryDihedral():
    system = [[FAI.get_atom(s='h'), FAI.get_atom(s='c',xyz=[1,1,3])],
            [FAI.get_atom(s='o',xyz=[0,1,1]), FAI.get_atom(s='s',xyz=[0,1,1])]]
    pt = [1,1,1]
    for i in system:
        for j in i:
            print(j.__dict__)
        print()
    idpts = [[0,1],[1,1]]
    fv = VaryDihedral(system,type='d',idpts=idpts,ratio=3.0)
    print(fv.v)
    print(fv.want,fv.type)
    print(fv.idpts)
    if not fv.nice:
        print(fv.info)
        return
    fv.run()
    print(fv.lower,fv.higher)
    print(fv.po, fv.pt, fv.pd)




class VaryZoom(VaryBond):
    """Zoom
    """
    def __init__(self,system,*args,**kwargs):
        if 'ipts' not in kwargs or kwargs['ipts'] is None:
            if 'idpts' not in kwargs or kwargs['idpts'] is None:
                idpts = []
                for i,mol in enumerate(system):
                    for j,at in enumerate(mol):
                        idpts.append([i,j])
                kwargs['idpts'] = idpts
  
        # be aware in here, we are calling its grandparent
        ParConfig.__init__(self,system,*args,**kwargs)
        if not self.nice: return
        if self.type != 'zoom':
            self.nice = False
            self.info += '\nFatal: configure error: not zooms variation'
            return
        getattr(self,'config_'+self.type)()




def test_class_VaryZoom():
    system = [[FAI.get_atom(s='h'), FAI.get_atom(s='c',xyz=[1,1,1])],
            [FAI.get_atom(s='o',xyz=[1,2,3]), FAI.get_atom(s='s',xyz=[0,1,1])]]

    for i in system:
        for j in i:
            print(j.__dict__)
        print()


    fv = VaryZoom(system,type='z',pd=[0,0,1],pts=[[0,0,1]])
    fv.run()
    print(fv.idpts)
    for i in range(5):
        print(fv.syspts[i])

    print(fv.num)
    for s in range(5):
        s = fv.sysnew[s]
        for m in s:
            for a in m:
                print(a.__dict__)
        print()
        print()

    xl = []
    yl = []
    zl = []
    for i in range(fv.num):
        if i % 30 == 0:
            for mol in fv.sysnew[i]:
                for at in mol:
                    xl.append(at.xyz[0])
                    yl.append(at.xyz[1])
                    zl.append(at.xyz[2])

    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(xl,yl,zl)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    plt.close()




class MolSampling:
    """sampling molecule in atomic spheres

    Args:
        system  : 2D n*nAtom : [ [Atom, ...],  ...]
        ipts    : uD [1,n]*[1,n]i  : [cnt, [cnt, radius], ...] counting starts at 0
        idpts   : 2D n*3     : [[mol-index, atom-index, radius], ...] start at 0

        mode (str|int) : all | partial | int | None, sample on setting, low priority

        seed (int): seed for random generation
        num (int) : number of generations, default 5

        sysnew  : 3D n*system
        axes (str) : axes to be sampled, x, y, z, default x*y*z*
            format:
                note: sign, * means either on + and -

                on axis both directions: x same as x*, xy same as x*y*, xyz same as x*y*z*
                otherwise: has to be clear, x+y-, y-, z*

        radius(float): uniformly all radius, highest priority
    """
    def __init__(self,system,*args,ipts=None,idpts=None,mode=None,num=None,seed=None,
                axes=None,radius=None,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system

        myipts = []
        myir = []
        if ipts is not None:
            if not isinstance(ipts,list):
                self.nice = False
                self.info = 'Fatal: wrong defined: not list type: ipts: {:}'.format(ipts)
                return
            for i in ipts:
                if isinstance(i,int):
                    myipts.append(i)
                    r = None
                elif isinstance(i,list) and len(i) == 2 and isinstance(i[0],int) \
                    and isinstance(i[1],(float,int)):
                        myipts.append(i[0])
                        r = i[1]
                else:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: ipts: {:}'.format(i)
                    return
                myir.append(r)
        
        myidpts = []
        myidr = []
        if idpts is not None:
            if not isinstance(idpts,list):
                self.nice = False
                self.info = 'Fatal: wrong defined: not list type: idpts: {:}'.format(idpts)
                return
            for i in idpts:
                bo = False
                if isinstance(i,list) and len(i) == 2:
                    if isinstance(i[0],int) and isinstance(i[1],int):
                        myidpts.append(i)
                        r = None
                    else:
                        bo = True
                elif isinstance(i,list) and len(i) == 3:
                    if isinstance(i[0],int) and isinstance(i[1],int) and isinstance(i[2],(float,int)):
                        myidpts.append(i[:2])
                        r = i[2]
                    else:
                        bo = True
                else:
                    bo = True
                if bo:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: idpts: {:}'.format(i)
                    return
                myidr.append(r)

        # get the help from ParConfig, set an additional temporary pts
        pc = ParConfig(system,pts=[[1,0,0]],ipts=myipts,idpts=myidpts)
        if not pc.nice:
            self.nice = False
            self.info = pc.info
            return

        self.ipts = pc.ipts

        self.idpts = []
        # note: idpts will combine with ipts, has to perform individual searching
        # priority: ipts > idpts
        if pc.idpts is not None and len(pc.idpts) != 0:
            myipts_id = [] if len(myipts) == 0 else pc.calc_pid(myipts)
            for cnt,j in enumerate(pc.idpts):
                if j in myipts_id:
                    rs = myir[myipts_id.index(j)]
                elif j in myidpts:
                    rs = myidr[myidpts.index(j)]
                else:
                    rs = None
                if rs is None: rs = self.system[j[0]][j[1]].r * 2
                self.idpts.append([j[0],j[1],rs])


        totlist = []
        for i,mol in enumerate(self.system):
            for j,at in enumerate(mol):
                totlist.append([i, j, at.r*2])

        # default mode
        if mode is None:
            if len(self.idpts) == 0: mode = 'all'

        if mode is not None and len(self.idpts) != 0:
            self.nice = False
            self.info = 'Fatal: conflict: mode and ipts & idpts cannot be set simultaneously'
            return
        
        # always set self.mode to None, detail in self.info
        self.mode = None
        if mode is None:
            pass
        elif isinstance(mode,str) and mode.lower() in ['all','full','a']:
            self.info = 'Note: mode all, sampling number: {:}'.format(len(totlist))
            self.idpts = totlist
        elif isinstance(mode,str) and mode.lower() in ['partial','random','p']:
            nm = np.random.randint(len(totlist))
            self.info = 'Note: mode partial, sampling number: {:}'.format(nm)
            refid = np.random.choice(range(len(totlist)),nm)
            self.idpts = [totlist[i] for i in refid]
        elif isinstance(mode,int):
            if mode > len(totlist):
                self.info = 'Warning: mode too large: reset to all'
                print(self.info)
                self.idpts = totlist
            else:
                refid = np.random.choice(range(len(totlist)),mode)
                self.idpts = [totlist[i] for i in refid]
        else:
            self.nice = False
            self.info = 'Fatal: wrong configure: mode: {:}'.format(mode)
            print(self.info)
            return

        if seed is None:
            np.random.seed()
            # extract a seed
            self.seed = np.random.get_state()[1][2]
        else:
            self.seed = seed
        np.random.seed(self.seed)

        self.num = 5 if num is None else num

        bo = False
        if axes is None:
            self.axes = 'x*y*z*'
        elif isinstance(axes,str):
            axes = axes.lower()
            if len(axes) == 1:
                if axes in 'xyz':
                    self.axes = axes + '*'
                else:
                    bo = True
            elif len(axes) == 2:
                if axes[0] in 'xyz' and axes[1] in 'xyz':
                    self.axes = axes[0] + '*' + axes[1] + '*'
                elif axes[0] in 'xyz' and axes[1] in '+-*':
                    self.axes = axes
                else:
                    bo = True
            elif len(axes) == 3:
                if axes[0] in 'xyz' and axes[1] in 'xyz' and axes[2] in 'xyz':
                    self.axes = axes[0] + '*' + axes[1] + '*' + axes[2] + '*'
                else:
                    bo = True
            elif len(axes) == 4:
                if axes[0] in 'xyz' and axes[1] in '+-*' \
                    and axes[2] in 'xyz' and axes[3] in '+-*':
                    if axes[0] == axes[2]:
                        bo = True
                    else:
                        self.axes = axes
                else:
                    bo = True
            elif len(axes) == 6:
                if axes[0] in 'xyz' and axes[1] in '+-*' \
                    and axes[2] in 'xyz' and axes[3] in '+-*' \
                        and axes[4] in 'xyz' and axes[5] in '+-*':
                    if axes[0] == axes[2] or axes[0] == axes[4] or axes[2] == axes[4]:
                        bo = True
                    else:
                        self.axes = axes
                else:
                    bo = True
            else:
                bo = True
        else:
            bo = True
        if bo:
            self.nice = False
            self.info = 'Fatal: wrong configure: axes: {:}'.format(axes)
            print(self.info)
            return
        
        if radius is None:
            pass
        elif isinstance(radius,(float,int)):
            self.radius = radius
            self.info = 'Note: uniformly set radius'
            for i in range(len(self.idpts)):
                self.idpts[i][2] = radius
        else:
            self.nice = False
            self.info = 'Fatal: wrong configure: radius: {:}'.format(radius)
            print(self.info)
            return



    def run(self):
        self.sysnew = []
        for i in range(self.num):
            mol = [[FAI.copy_atom(i) for i in j] for j in self.system]
            self.sysnew.append(mol)
        
        # in-place update
        if 'x+' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[0] += np.random.rand() * p[2]
        elif 'x-' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[0] -= np.random.rand() * p[2]
        elif 'x*' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    tmp = np.random.choice([-1,0,1],1)[0]
                    mol[p[0]][p[1]].xyz[0] += np.random.rand() * p[2] * tmp
        
        if 'y+' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[1] += np.random.rand() * p[2]
        elif 'y-' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[1] -= np.random.rand() * p[2]
        elif 'y*' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    tmp = np.random.choice([-1,0,1],1)[0]
                    mol[p[0]][p[1]].xyz[1] += np.random.rand() * p[2] * tmp
        
        if 'z+' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[2] += np.random.rand() * p[2]
        elif 'z-' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    mol[p[0]][p[1]].xyz[2] -= np.random.rand() * p[2]
        elif 'z*' in self.axes:
            for mol in self.sysnew:
                for p in self.idpts:
                    tmp = np.random.choice([-1,0,1],1)[0]
                    mol[p[0]][p[1]].xyz[2] += np.random.rand() * p[2] * tmp




def test_class_MolSampling():
    system = [[FAI.get_atom(s='h'), FAI.get_atom(s='c',xyz=[1,1,1])],
            [FAI.get_atom(s='o',xyz=[1,2,3]), FAI.get_atom(s='s',xyz=[0,1,1])]]

    for i in system:
        for j in i:
            print(j.__dict__)
        print()

    fv = MolSampling(system)
    assert fv.mode is None
    assert 'all' in fv.info

    fv = MolSampling(system,mode='p')
    assert fv.mode is None
    assert 'partial' in fv.info

    fv = MolSampling(system,mode=2)
    assert fv.mode is None
    assert len(fv.idpts) == 2

    fv = MolSampling(system,mode=1000)
    assert fv.mode is None
    assert len(fv.idpts) == 4
    assert fv.axes == 'x*y*z*'
    
    fv = MolSampling(system,mode=1000,ipts=[0,1])
    assert not fv.nice

    fv = MolSampling(system,idpts=[[0,0,1]],axes='y-x+z*',radius=2)
    assert fv.nice
    assert np.allclose(fv.idpts, [[0,0,2]])
    fv.run()
    print('printing')
    for s in range(5):
        for m in fv.sysnew[s]:
            for a in m:
                print(a.__dict__)
        print()
        print()


    fv = MolSampling(system,ipts=[[0,2],2],idpts=[[0,0,1]])
    assert fv.nice
    assert len(fv.idpts) == 2
    print(fv.idpts)
    fv.run()

    print(fv.mode)
    print(fv.seed)
    print(fv.idpts)
    print(fv.info)
    print(fv.num)

    for s in range(5):
        for m in fv.sysnew[s]:
            for a in m:
                print(a.__dict__)
        print()
        print()




class FixBonds:
    """
    Args:
        system  : 2D n*nAtom : [ [Atom, ...],  ...]
        sysnew  : 3D n*system

                system is used for reference to calculate bcon, low priority
                if it is omitted, sysnew[0] will be used instead

        bcon    : 2D [[atom-i, atom-j, dt], ...] : List[[int,int,float], ...]
                  dt is the correctness, can be omitted
                  format on user inputs, starts at 1
                  attributes on codes use, starts at 0

                note: the index is counted as the molecule as a whole,
                e.g., for system [[a,b,c],[s,t]], their indexes will be taken
                as an integral, [[1,2,3],[4,5]], thus, if bonded info as,
                [a-b, a-c, a-t, s-t-0.8], then bcon should be set to:

                *) [[1,2], [1,3], [1,5], [4,5,0.8]], user inputs, starts at 1
                **) [[0,1], [0,2], [0,4], [3,4,0.8]], codes use, starts at 0
                ***) [[0,0,0,1], [0,0,0,2], [0,0,1,0], [1,0,1,1,0.8]], hard-coded

                take the last one, [1,0,1,1,0.8]], for explanation
                i)   1  (1/5, or 1/4) => mol-index for atom s
                ii)  0  (2/5, or 2/4) => atom-index for atom s
                iii) 1  (3/5, or 3/4) => mol-index for atom t
                iv)  1  (4/5, or 4/4) => atom-index for atom t
                v)   0.8(5/5), correctness, can be omitted

                caution: codes use and hard-coded are not invisible,
                instead, user inputs can be added as the supplementary
                please refer to bmode and set it correctly to 'add'.

                besides, filters can be set to remove entries of bcon, only the
                correct settings will be taken, no worries about input errors

        bmode (str|None): add,added,addition,a | part,separate,p,s | None, lowest priority
                for option add, which has been discussed in bcon.

                for part, which means the bond perception is only performed on individual
                molecules, thus if bcon is not None and part is turned on, warning info
                will be printed for wrong inputs and they will be thrown away,
                correct inputs will be added instead

        filters : same format as bcon, implicitly filter list, highest priority
                  format on user inputs, starts at 1

        sysnew will be updated filtered results after self.run() is executed

    Reference:
        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Equation:
        0.8 <= dij <= ri + rj + 0.4 + correctness
    """
    def __init__(self,system,*args,sysnew=None,bcon=None,bmode=None,filters=None,
                userinputs=None,**kwargs):
        self.nice = True
        self.info = ''
        if sysnew is None:
            self.sysnew = system
            self.system = system[0]
        else:
            self.system = system
            self.sysnew = sysnew

        self.userinputs = False if userinputs is False else True

        if bcon is None:
            self.bcon = self.calc_bond_perception(self.system)
        elif userinputs:
            self.bcon = self.check_user_inputs(bcon)
            if not self.nice: return
        else:
            self.bcon = bcon

        if bmode is None:
            self.bmode = None
        elif bmode.lower() in ['add','added','addition','a']:
            self.bmode = 'add'
            if bcon is not None:
                auto = self.calc_bond_perception(self.system)
                for con in auto:
                    bo = False
                    for ndx in self.bcon:
                        if con[0] == ndx[0] and con[1] == ndx[1]:
                            bo = True
                            break
                    if not bo:
                        self.bcon.append(con)
        elif bmode.lower() in ['part','separate','p','s']:
            self.bmode = 'part'
            chk = []
            if bcon is not None:
                tmplist = [len(i) for i in self.system]
                acclist = [0,]
                tot = 0
                for i in tmplist:
                    tot += i
                    acclist.append(tot)
                for con in self.bcon:
                    i = 1
                    while i < len(acclist):
                        if con[0] < acclist[i]:
                            break
                        i += 1
                    if con[1] >= acclist[i-1] and con[1] < acclist[i]:
                        chk.append(con)
                    else:
                        print('Warning: bmode part filters: removing {:}'.format(con))
            self.bcon = self.calc_bond_perception(self.system,boall=False)
            for i in chk: self.bcon.append(i)
        else:
            print('Warning: wrong defined: reset to None: bmode: {:}'.format(bmode))
            self.bmode = None

        if filters is None:
            self.filters = None
        else:
            # care, user inputs, starts at 1
            # process user inputs, collect valid filters
            self.filters = []
            for con in filters:
                bo = False
                for ndx in self.bcon:
                    if con[0]-1 == ndx[0] and con[1]-1 == ndx[1]:
                        bo = True
                        break
                if bo:
                    self.bcon.remove(ndx)
                    self.filters.append(con)
    


    def check_user_inputs(self,bcon):
        rstlist = []
        atmax = sum([len(i) for i in self.system])
        for ndx in bcon:
            if len(ndx) == 2 or len(ndx) == 3:
                if max(ndx[:2]) > atmax or min(ndx[:2]) <= 0:
                    self.nice = False
                    self.info = 'Warning: wrong defined: idpts: {:}'.format(ndx)
                    return
                if len(ndx) == 2:
                    tmp = [ndx[0]-1, ndx[1]-1]
                else:
                    tmp = [ndx[0]-1, ndx[1]-1, ndx[2]]
                rstlist.append(tmp)
            else:
                self.nice = False
                self.info = 'Warning: wrong format : reset to all: {:}'.format(ndx)
                return
        return rstlist



    def calc_bond_perception(self,system,boall=True):
        if len(system) == 1 and len(system[0]) == 1: return []
        bcon = []
        if boall:
            # squeeze system to 1D
            systmp = [at for mol in system for at in mol]
            for i,ai in enumerate(systmp[:-1]):
                k = i + 1
                for j,aj in enumerate(systmp[k:]):
                    dx = aj.xyz[0] - ai.xyz[0]
                    dy = aj.xyz[1] - ai.xyz[1]
                    dz = aj.xyz[2] - ai.xyz[2]
                    dd = dx*dx + dy*dy + dz*dz
                    t = ai.r + aj.r + 0.4
                    if dd > 0.8 and dd < t*t:
                        bcon.append([i,j+k])
        else:
            tot = 0
            for systmp in system:
                for i,ai in enumerate(systmp[:-1]):
                    k = i + 1
                    for j,aj in enumerate(systmp[k:]):
                        dx = aj.xyz[0] - ai.xyz[0]
                        dy = aj.xyz[1] - ai.xyz[1]
                        dz = aj.xyz[2] - ai.xyz[2]
                        dd = dx*dx + dy*dy + dz*dz
                        t = ai.r + aj.r + 0.4
                        if dd > 0.8 and dd < t*t:
                            bcon.append([i+tot,j+k+tot])
                tot += len(systmp)
        return bcon



    def calc_reflist(self,system,bcon):
        # split total index to individual index
        # corresponding to hard-coded info
        tmplist = [len(i) for i in system]
        acclist = [0,]
        tot = 0
        for i in tmplist:
            tot += i
            acclist.append(tot)
        reflist = []
        for ndx in bcon:
            for i,nm in enumerate(acclist):
                if ndx[0] < nm:
                    mi = i - 1
                    ai = ndx[0] - acclist[i-1]
                    break
            for i,nm in enumerate(acclist):
                if ndx[1] < nm:
                    mj = i - 1
                    aj = ndx[1] - acclist[i-1]
                    break
            if len(ndx) == 2:
                reflist.append([mi,ai,mj,aj,0.0])
            else:
                reflist.append([mi,ai,mj,aj,ndx[2]])
        return reflist



    def run(self):
        if len(self.bcon) == 0:
            print('Warning: no filtration: bcon is empty')
            return

        self.reflist = self.calc_reflist(self.system,self.bcon)

        sysrst = []
        for system in self.sysnew:
            bo = True
            for ndx in self.reflist:
                ai = system[ndx[0]][ndx[1]]
                aj = system[ndx[2]][ndx[3]]
                dx = aj.xyz[0] - ai.xyz[0]
                dy = aj.xyz[1] - ai.xyz[1]
                dz = aj.xyz[2] - ai.xyz[2]
                dd = dx*dx + dy*dy + dz*dz
                t = ai.r + aj.r + 0.4 + ndx[4]
                if dd < 0.64 or dd > t*t:
                    bo = False
                    break
            if bo:
                sysrst.append(system)
        self.num = len(sysrst)
        self.sysnew = sysrst



def test_class_FixBonds():
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000
        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325

        H           -1.081          -1.003          1.416
        H            0.576           -0.336         1.933
        C           -0.329            1.632         2.341
        O            1.268            0.242         1.064
        N           -0.017            0.349         4.373
        10           1.331            2.040         3.305
        H            0.155           -0.527         2.467
        H            1.118            0.107         3.318
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    rf = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    rf.run()

    system = FAI.guess_atoms_for_system(rf.system)

    fv = MolSampling(system,mode='all',axes='y',num=1000)
    if not fv.nice:
        print(fv.info)
        return
    fv.run()
    print(fv.num)

    bcon = []
    fb = FixBonds(fv.system,sysnew=fv.sysnew,bmode='p')
    if not fb.nice:
        print(fb.info)
        exit()
    print(fb.bcon,fb.bmode,fb.filters)
    fb.run()

    print(len(fb.sysnew))
    for i in range(min(len(fb.sysnew),2)):
        for mol in fb.sysnew[i]:
            for at in mol:
                print(at.__dict__)
        print()
        print()




class FixNonBonds(FixBonds):
    """
    Args:
        nbcon : 2D n*2i : nonbonded index list, performed filtration on

        sysnew will be updated filtered results after self.run() is executed
    """
    def __init__(self,system,*args,covalent=None,**kwargs):
        super().__init__(system,*args,**kwargs)
        if not self.nice: return
        nbcon = kwargs['nbcon'] if 'nbcon' in kwargs else None
        if nbcon is None:
            tot = sum([len(i) for i in self.system])
            self.nbcon = []
            for i in range(tot-1):
                for j in range(i+1,tot):
                    bo = False
                    for ndx in self.bcon:
                        if ndx[0] == i and ndx[1] == j:
                            bo = True
                            break
                    if not bo:
                        self.nbcon.append([i,j])
        elif self.userinputs:
            self.nbcon = self.check_user_inputs(nbcon)
            if not self.nice: return
        else:
            self.nbcon = nbcon

        self.covalent = covalent



    def run(self):
        if len(self.nbcon) == 0:
            print('Warning: no filtration: nbcon is empty')
            return

        self.reflist = self.calc_reflist(self.system,self.nbcon)

        sysrst = []
        hh = 0.74*0.74 if self.covalent is None else self.covalent*self.covalent
        for system in self.sysnew:
            bo = True
            for ndx in self.reflist:
                ai = system[ndx[0]][ndx[1]]
                aj = system[ndx[2]][ndx[3]]
                dx = aj.xyz[0] - ai.xyz[0]
                dy = aj.xyz[1] - ai.xyz[1]
                dz = aj.xyz[2] - ai.xyz[2]
                dd = dx*dx + dy*dy + dz*dz
                if dd < hh:
                    bo = False
                    break
            if bo:
                sysrst.append(system)
        self.num = len(sysrst)
        self.sysnew = sysrst




def test_class_FixNonbonds():
    ftxt = """
        C     0.000     0.000       0.000
        C     0.000     0.000       1.500
        O     1.226     0.000      -0.490
        N    -1.370    -0.003       2.109
        O    -1.934     1.053       2.350
        O    -1.912    -1.060       2.391
        H     0.595     0.880       1.848
        H     0.601    -0.876       1.849
        H    -0.610     0.891      -0.335
        H    -0.611    -0.891      -0.335
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    rf = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    rf.run()

    system = FAI.guess_atoms_for_system(rf.system)

    fv = MolSampling(system,ipts=[0,6],axes='x-',num=1000,radius=4)
    if not fv.nice:
        print(fv.info)
        return
    print(fv.idpts)
    fv.run()
    print(fv.num)
    fnb = FixNonBonds(system,bcon=[[1,7]],sysnew=fv.sysnew)
    assert fnb.nice
    print('bcon =',fnb.bcon)
    print('nbcon=',fnb.nbcon)
    fnb.run()
    print('second ',len(fnb.sysnew))




class GenBonds(VaryBond,FixNonBonds):
    def __init__(self,system,*args,**kwargs):
        bo = False
        if 'idpts' not in kwargs or kwargs['idpts'] is None:
            bo = True
            kwargs['idpts'] = [[0,0]]
        VaryBond.__init__(self,system,**kwargs)
        if not self.nice:
            print(self.info)
            return
        if bo:
            self.idpts = None
            self.calc_pars()

        # be aware in here, system is added an additional dimension
        FixNonBonds.__init__(self,[system],**self.kwargs)
        if not self.nice:
            print(self.info)
            return



    def calc_pars(self):
        self._calc_accumulates()

        # format:
        #       dual|single     num-pts     [kwargs, kwargs ...]
        n = sum([len(i) for i in self.system])
        self.totpars = {'dual':{}, 'single':{}}
        self.totnews = {'dual':{}, 'single':{}}
        for p in range(1,n-1):
            self.totpars['dual'][p] = []
            self.totpars['single'][p] = []

            self.totnews['dual'][p] = []
            self.totnews['single'][p] = []

        if self.idpts is None:
            bopo = True if 'po' in self.kwargs and self.kwargs['po'] is not None else False
            bopd = True if 'pd' in self.kwargs and self.kwargs['pd'] is not None else False
            if (not bopo) and (not bopd):
                for i in range(n):
                    for j in range(n):
                        if j == i: continue

                        ai = self.calc_indexs(i)
                        aj = self.calc_indexs(j)

                        # only direction point
                        # be aware, idpts should be in 2D
                        self.totpars['dual'][1].append(
                            {
                                'idpo': ai,
                                'idpd': aj,
                                'idpts': [aj],
                                'inc': 0.03,
                                'ratio': 3,
                            }
                        )
                        self.totpars['dual'][1].append(
                            {
                                'idpo': ai,
                                'idpd': aj,
                                'idpts': [aj],
                                'inc': 0.03,
                                'ratio': 0.3,
                            }
                    )

                        s = [k for k in range(n) if k != i and k != j]
                        for r in range(1,n-1):
                            for ndx in itertools.combinations(s,r):
                                mm = self.calc_indexs(ndx)
                                if r != n - 2:
                                    # Care! shallow copy of aj
                                    mm_with_aj = [aj, ]
                                    mm_with_aj.extend(mm)

                                    # dual when contains direction point
                                    self.totpars['dual'][r+1].append(
                                        {
                                            'idpo': ai,
                                            'idpd': aj,
                                            'idpts': mm_with_aj,
                                            'inc': 0.03,
                                            'ratio': 3,
                                        }
                                    )
                                    self.totpars['dual'][r+1].append(
                                        {
                                            'idpo': ai,
                                            'idpd': aj,
                                            'idpts': mm_with_aj,
                                            'inc': 0.03,
                                            'ratio': 0.3,
                                        }
                                    )

                                # single when no direction point
                                self.totpars['single'][r].append(
                                    {
                                        'idpo': ai,
                                        'idpd': aj,
                                        'idpts': mm,
                                        'inc': 0.03,
                                        'ratio': 3,
                                    }
                                )



    def run(self):
        pars = self.totpars['dual']
        for nm in pars:
            for p in pars[nm]:
                VaryBond.__init__(self,self.system,**p)
                VaryBond.run(self)

                if self.bcon is not None and len(self.bcon) != 0:
                    FixBonds.run(self)
                if self.nbcon is not None and len(self.nbcon) != 0:
                    FixNonBonds.run(self)

                self.totnews['dual'][nm].append(self.sysnew)



        pars = self.totpars['single']
        for nm in pars:
            for p in pars[nm]:
                VaryBond.__init__(self,self.system,**p)
                VaryBond.run(self)

                if self.bcon is not None and len(self.bcon) != 0:
                    FixBonds.run(self)
                if self.nbcon is not None and len(self.nbcon) != 0:
                    FixNonBonds.run(self)

                self.totnews['single'][nm].append(self.sysnew)



    def calc_indexs(self,atoms):
        """calculate individual index to [mol-index, atom-index]

        Args:
            atoms: int | 1D 1*ni

        Return:
            1D 1*2i     [int, int]          for int
            2D n*2i     [[int,int], ...]    Otherwise
        """
        if isinstance(atoms,int):
            for i,nm in enumerate(self._acclist):
                if atoms < nm:
                    mi = i - 1
                    ai = atoms - self._acclist[i-1]
                    break
            return [mi,ai]

        reflist = []
        for atnum in atoms:
            for i,nm in enumerate(self._acclist):
                if atnum < nm:
                    mi = i - 1
                    ai = atnum - self._acclist[i-1]
                    break
            reflist.append([mi,ai])
        return reflist


    def _calc_accumulates(self):
        tmplist = [len(i) for i in self.system]
        self._acclist = [0,]
        tot = 0
        for i in tmplist:
            tot += i
            self._acclist.append(tot)




def test_class_GenBonds():
    ai = AtomInfo()
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000

        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()

    system = ai.guess_atoms_for_system(RF.system)

    GB = GenBonds(system)
    GB.run()
    if not GB.nice:
        print(GB.info)
        return

    pars = GB.totpars['dual']
    news = GB.totnews['dual']
    for nm in pars:
        kws = pars[nm]
        totmols = news[nm]

        ndx = 1
        print(kws[ndx])
        mols = totmols[ndx]
        print(len(mols))
        for j in range(min(len(mols),2)):
            print('num -->',j+1)
            for mol in mols[j]:
                for at in mol:
                    print(at.s,at.xyz)
                print()
            print('\n')
        print('\n\n')




class SaveFileSystemSingle:
    """opposite operation to ReadFile

    Args:
        system : 3D System[ Mol[Atom, ...], Mol[Atom, ...] , ... ]
        ftype  : Output file type
            format
                txt file
                xsf file
                xyz file
        fname : file to be saved, warning, overwritten may happen
    """
    def __init__(self,system,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system
        if len(self.system) == 0:
            self.nice = False
            self.info = 'Warning: no inputs: system'
            return


        # Rule
        # self.ftype always takes the precedence
        # a) if self.fname has the period
        #       1) if its file extension matchs with self.ftype
        #          everything is fine
        #       2) if not, append real ftype on it
        # b) append real ftype
        self.ftype = None
        if 'ftype' in kwargs and kwargs['ftype'] is not None:
            self.ftype = kwargs['ftype'].lower()

        self.fname = None
        if 'fname' in kwargs and kwargs['fname'] is not None:
            if len(kwargs['fname'].split()) != 0:
                self.fname = kwargs['fname']

        if self.fname is not None and self.ftype is None:
            # guess ftype from fname
            ndx = self.fname.rfind('.')
            if ndx != -1:
                ext = self.fname[ndx:]
                if ext == '.':
                    self.fname = self.fname[:-1]
                elif ext == '.txt':
                    self.ftype = 'txt'
                elif ext == '.xsf':
                    self.ftype = 'xsf'
                elif ext == '.xyz':
                    self.ftype = 'xyz'
                elif ext == '.com':
                    self.ftype = 'com'

        if self.ftype is None: self.ftype = 'txt'
        if self.ftype not in ['txt','xsf','xyz','com']:
            self.nice = False
            print('Warning: not support < {:} >'.format(self.ftype))
            self.info = 'Warning: not support < {:} >'.format(self.ftype)
            return 

        if self.fname is None: self.fname = 'system.txt'
        # keep dot conversion in original fname, if it has
        ndx = self.fname.rfind('.')
        if ndx == -1:
            self.fname = self.fname + '.' + self.ftype
        else:
            ext = self.fname[ndx:]
            if ext != '.' + self.ftype:
                self.fname = self.fname + '.' + self.ftype


        # alias final chosen function
        self.fn = getattr(self,'yield_'+self.ftype)



    def run(self):
        with open(self.fname,'wt') as f:
            f.write(self.fn(self.system))



    def yield_xsf(self,system):
        fout = '#\n\nATOMS\n'
        for mol in system:
            for at in mol:
                fout += '{:2} {:>15.8f} {:>15.8f} {:>15.8f}   1.0  1.0  1.0\n'.format(at.s,*at.xyz)
        fout += '\n\n'
        return fout



    def yield_txt(self,system):
        fout = ''
        for mol in system:
            for at in mol:
                fout += '{:2} {:>15.8f} {:>15.8f} {:>15.8f}\n'.format(at.s,*at.xyz)
        fout += '\n\n'
        return fout



    def yield_xyz(self,system):
        fout = '{:}\n'.format(sum([len(i) for i in system]))
        fout += 'Properties=species:S:1:pos:R:3 energy=0.0\n'
        for mol in system:
            for at in mol:
                fout += '{:2} {:>15.8f} {:>15.8f} {:>15.8f}\n'.format(at.s,*at.xyz)
        fout += '\n\n'
        return fout
    


    def yield_com(self,system):
        fout = '#\n\nTitle\n\n0 1\n'
        for mol in system:
            for at in mol:
                fout += '{:2} {:>15.8f} {:>15.8f} {:>15.8f}\n'.format(at.s,*at.xyz)
        fout += '\n\n'
        return fout




def test_class_SaveFileSystemSingle():
    ai = AtomInfo()
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000

        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()
    
    system = ai.guess_atoms_for_system(RF.system)

    SFS = SaveFileSystemSingle(system)
    if not SFS.nice:
        print(SFS.info)
        return

    print(SFS.fname, SFS.ftype)
    SFS.run()




class SaveFileSystemMany(SaveFileSystemSingle):
    def __init__(self,system,*args,**kwargs):
        super().__init__(system,*args,**kwargs)

        self.showcase = False
        if 'showcase' in kwargs and kwargs['showcase'] is True:
            self.showcase = True



    def run(self):
        if self.showcase:
            self.run_all()
            return
        with open(self.fname,'wt') as f:
            for each in self.system:
                f.write(self.fn(each))



    def run_all(self):
        name = self.fname[:self.fname.rfind('.')] + '.com'
        with open(name,'wt') as f:
            f.write('#\n\nMOSS showcase\n\n0 1\n')
            for each in self.system:
                for mol in each:
                    for at in mol:
                        l = '{:2} {:>15.8f} {:>15.8f} {:>15.8f}\n'.format(at.s,*at.xyz)
                        f.write(l)
            f.write('\n\n\n')




def test_class_SaveFileSystemMany():
    ai = AtomInfo()
    ftxt = """
        C           -0.225          -0.346          1.167
        1           -0.917           0.728          1.999
        2            0.000           0.000          0.000

        N           -0.978          -0.007          3.141
        O           -1.950          -0.786          3.325
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()
    
    system = ai.guess_atoms_for_system(RF.system)


    SFM = SaveFileSystemMany([system for i in range(5)])
    if not SFM.nice:
        print(SFM.info)
        return

    print(SFM.fname, SFM.ftype)
    SFM.run()




def save_genbonds():
    ai = AtomInfo()
    ftxt = """
        C    -2.37478119    3.11200745    0.06229865
        O    -1.75709606    3.85995114    1.11298667
        F    -3.71663176    3.22455427    0.15858084
        Cl   -1.91671366    1.41827100    0.20025150
        N    -1.93120844    3.63524554   -1.23784177
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()
    
    system = ai.guess_atoms_for_system(RF.system)


    GB = GenBonds(system)
    if not GB.nice:
        print(GB.info)
        return
    GB.run()


    pars = GB.totpars['dual']
    news = GB.totnews['dual']
    for nm in pars:
        for cnt,kw in enumerate(pars[nm]):
            systems = news[nm][cnt]

            fname = 'sample-dual-{:}-{:}-{:}'.format(nm,kw['idpo'][1]+1,kw['idpd'][1]+1)
            SFS = SaveFileSystemMany(systems,fname=fname,showcase=True)
            SFS.run()


    pars = GB.totpars['single']
    news = GB.totnews['single']
    for nm in pars:
        for cnt,kw in enumerate(pars[nm]):
            systems = news[nm][cnt]

            fname = 'sample-single-{:}-{:}-{:}'.format(nm,kw['idpo'][1]+1,kw['idpd'][1]+1)
            SFS = SaveFileSystemMany(systems,fname=fname,showcase=True)
            SFS.run()


def gen_bonds_overall():
    ai = AtomInfo()
    ftxt = """
        # Henry reference
        C     0.000     0.000       0.000
        C     0.000     0.000       1.500
        O     1.226     0.000      -0.490
        N    -1.370    -0.003       2.109
        O    -1.934     1.053       2.350
        O    -1.912    -1.060       2.391
        H     0.595     0.880       1.848
        H     0.601    -0.876       1.849
        H    -0.610     0.891      -0.335
        H    -0.611    -0.891      -0.335
    """
    fp = tempfile.TemporaryFile()
    fp.write(ftxt.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    RF = ReadFileMultiple(fp.name)
    # fp will be destoried inside open!
    RF.run()
    
    system = ai.guess_atoms_for_system(RF.system)


    GB = GenBonds(system)
    if not GB.nice:
        print(GB.info)
        return
    GB.run()


    overall = []

    pars = GB.totpars['dual']
    news = GB.totnews['dual']
    for nm in pars:
        for cnt,kw in enumerate(pars[nm]):
            systems = news[nm][cnt]
            overall.extend(systems)

    pars = GB.totpars['single']
    news = GB.totnews['single']
    for nm in pars:
        for cnt,kw in enumerate(pars[nm]):
            systems = news[nm][cnt]
            overall.extend(systems)

    print('number of molecules ==>',len(overall))
    fname = 'sample-overall'
    SFS = SaveFileSystemMany(overall,fname=fname)
    SFS.run()


gen_bonds_overall()

DONE = 'all done'
print(DONE)