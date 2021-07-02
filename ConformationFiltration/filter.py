#!/usr/bin/env python3

import os
import sys
import math
import argparse
import matplotlib.pyplot as plt
import tempfile


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
    'version 2.60 : add FILEFORMAT explanation',
    'version 2.70 : fix CrossFiltration reflist',
    'version 2.80 : fix ReadFile xsf',
    'version 2.90 : make ReadFile xsf more flexible',
    'version 3.00 : add save to xyz file',
    'version 3.10 : make fname and ftype more compatible in save to file',
    'version 3.20 : add read xyz file',
    'version 3.21 : update FILEFORMAT for xyz',
    'version 3.30 : make ReadFile more compatible',
    'version 3.3.0  : BC, add class AtomInfo',
    'version 3.3.1  : refine ReadFile',
    'version 3.3.2  : refine SaveFile',
    'version 3.3.3  : refine BondPerception',
    'version 3.3.4  : add AnglePerception',
    'version 3.3.5  : class perceptions are done',
    'version 3.3.6  : refine Filtration',
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
"""


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
        self.atominfo = []
        for line in self.ATOM_PROPERTY.split('\n')[1:]:
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
            self.atominfo.append([sign,number,radius,mass,name])


FAI = AtomInfo()


class ReadFile:
    """
    Args:
        file (str): input file name
        ext (str): txt | xsf | xyz
        debug (bool): whether printout more info

    Attributes:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        energy : 1D List[float]  :   None means not exist
    """
    def __init__(self,file,ext=None,debug=True,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.file = file

        # decide file format
        if ext is None:
            ndx = file.rfind('.')
            if ndx == -1 or ndx + 1 >= len(file):
                self.ext = 'txt'
            else:
                self.ext = file[ndx+1:].lower()
        else:
            self.ext = ext
        
        if self.ext not in ['txt','xsf','xyz']:
            self.nice = False
            self.info = 'Fatal: file format not support: {:}'.format(ext)
            return
        
        self.debug = True if debug is True else False



    def run(self):
        """process input file

        require:
            atomtype entry can be mixed, which means,

            a) 1  x  y  z
            b) H  x  y  z

            they are equivalent
        """
        self.system = []
        self.energy = []

        prolist,enelist,errlist = getattr(self,'read_'+self.ext)()
        if len(prolist) <= 1: return

        if self.debug:
            for i in errlist:
                print('Warning: ignoring: {:}: {:}'.format(i[0],i[1]))
        
        self.system.append(prolist[0])
        self.energy.append(enelist[0])

        # format: 2D str: [ [sign, number, name], ... ]
        atominfo = [[i[0], str(i[1]), i[4]] for i in FAI.atominfo]

        ndxlist = []
        for i in prolist[0]:
            atype = i[0].capitalize()
            bo = True
            for ndx in atominfo:
                if atype in ndx:
                    bo = False
                    atype = ndx[0]
                    break
            if bo: atype = i[0]
            ndxlist.append(atype)

        nats = len(ndxlist)
        errlist = []
        for n,mol in enumerate(prolist[1:]):
            if len(mol) != nats:
                info = 'Warning: ignoring: error: number of atoms'
                errlist.append([info,mol])
                continue

            # mixed atomtype may occur
            # always make sure mol be in real atomtype
            # be aware of how index in mol is changed, python pointer
            bo = True
            for cnt,atom in enumerate(mol):
                if atom[0] != ndxlist[cnt]:
                    bo = False
                    atype = atom[0].capitalize()
                    for ndx in atominfo:
                        if atype in ndx:
                            bo = True
                            atype = ndx[0]
                            break
                    if bo:
                        if atype == ndxlist[cnt]:
                            atom[0] = atype
                        else:
                            bo = False
                    if not bo:
                        break
            if bo:
                # be careful with n starts
                self.energy.append(enelist[n+1])
                self.system.append(mol)
            else:
                info = 'Warning: ignoring: error: not cooresponded'
                errlist.append([info,mol])
    
        if self.debug:
            for i in errlist:
                print(i[0])
                for j in i[1]: print(j)
                print()



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
        
        prolist = []
        enelist = []
        errlist = []
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
                ene = None
                ltmp = mol[0][0].replace('=',' ').split()
                if len(ltmp) >= 2:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass

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
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist



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

        prolist = []
        enelist = []
        errlist = []
        for mol in promol:
            # check whether energy exist or not
            ene = None
            if mol[0][0][0] == '#':
                ltmp = mol[0][0].replace('=',' ').split()
                if len(ltmp) >= 2:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass
                mol = mol[1:]

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
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist



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

        prolist = []
        enelist = []
        errlist = []
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
                ene = None
                # check whether energy exist or not
                ltmp = mol[1][0].replace('=',' ').split()
                if len(ltmp) > 1:
                    try:
                        ene = float(ltmp[-1])
                    except ValueError:
                        pass

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
                errlist.append([errnum,errline])
            else:
                prolist.append(ls)
                enelist.append(ene)

        return prolist, enelist, errlist




class SaveFile:
    """opposite operation to ReadFile

    Args:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        energy : 1D List[float]  :   None means not exist

        ftype (str): Output file type: txt | xsf | xyz
        fname (str): file to be saved, warning, overwritten may happen
    
    Note:
        1) number of atoms in system will not be crossly checked, which means
        for: system[ mol(8),  mol(5),  mol(20), ...],
        it can be successfully saved

        2) saved atomtype will always be real atomtype, if it is valid
    """
    def __init__(self,system,energy=None,ftype=None,fname=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system if isinstance(system,list) else []
        self.energy = energy if energy is not None and isinstance(energy,list) else []

        if len(self.system) == 0:
            self.nice = False
            self.info = 'Fatal: no inputs'
            return

        # Rule
        # self.ftype always takes the precedence
        # a) if self.fname has the period
        #       1) if its file extension matchs with self.ftype
        #          everything is fine
        #       2) if not, append real ftype on it
        # b) append real ftype
        self.ftype = None if ftype is None else ftype.lower()

        self.fname = None
        if fname is not None and isinstance(fname,str) and len(fname.split()) != 0:
            self.fname = fname.strip()

        if self.fname is not None and self.ftype is None:
            # guess ftype from fname
            ndx = self.fname.rfind('.')
            if ndx != -1:
                ext = self.fname[ndx:]
                if ext == '.':
                    self.fname = self.fname[:-1]
                else:
                    self.ftype = ext.replace('.','')

        if self.ftype is None: self.ftype = 'txt'
        if self.ftype not in ['txt','xsf','xyz']:
            self.nice = False
            self.info = 'Fatal: not support: {:}'.format(self.ftype)
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

        # get real atomtype list
        # format: 2D: [ [sign, number-int, number-str, name], ... ]
        atominfo = [[i[0], i[1], str(i[1]), i[4]] for i in FAI.atominfo]
        self.atypelist = []
        for mol in self.system:
            ls = []
            for at in mol:
                atype = at[0].capitalize() if isinstance(at,str) else at[0]
                for ndx in atominfo:
                    if atype in ndx:
                        atype = ndx[0]
                        break
                ls.append(atype)
            self.atypelist.append(ls)



    def run(self):
        fout = getattr(self,'save_'+self.ftype)()
        with open(self.fname,'wt') as f: f.write(fout)



    def save_xsf(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                fout += '# {:}\n\nATOMS\n'.format(self.energy[ndx])
            else:
                fout += '#\n\nATOMS\n'
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:3} {:>10} {:>10} {:>10}   1.0  1.0  1.0\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout



    def save_txt(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                fout += '#  {:}\n'.format(self.energy[ndx])
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:<2} {:>15} {:>15} {:>15}\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout



    def save_xyz(self):
        fout = ''
        for ndx,mol in enumerate(self.system):
            fout += '{:}\n'.format(len(mol))
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                fout += 'Properties=species:S:1:pos:R:3 energy={:}\n'.format(self.energy[ndx])
            else:
                fout += 'Properties=species:S:1:pos:R:3 energy=0.0\n'
            for cnt,at in enumerate(mol):
                atype = self.atypelist[ndx][cnt]
                fout += '{:2} {:>12} {:>12} {:>12}\n'.format(atype,*at[1:])
            fout += '\n\n'
        return fout




def test_class_ReadFile_and_SaveFile():
    txtfile = """
        # -223997.547055016209755
        6 -0.13274  0.09930  1.31039
        8  0.00000  0.00000  0.00000
        7  0.42454  1.52049  3.33357
        1 -0.48212  2.28202  1.54152

        C -0.13274  0.09930  1.31039
        8  0.00000  0.00000  0.00000
        7  0.42454  1.52049  3.33357    2.333   good
        H  0.49216 -0.64485  1.82439

        # -223996.527301928327235
        carbon      -0.09295  0.02907  1.30111  still good
        Oxygen      -0.01327 -0.06802 -0.01374
        7            0.64812  1.44918  3.29189
        Hydrogen    -0.45954  2.22518  1.62322


        6 -0.02976  0.02692  1.34829    ignored, due to number of atoms
        6  0.37245  1.46435  1.85092
        8  0.04793 -0.07239  0.03349
        7  0.70929  1.43837  3.32682
        1 -0.40465  2.19862  1.65751


        6 -0.12917  0.11157  1.31695    ignored, due to number of atoms
        8  0.00308  0.00995  0.00669


        6 -0.12917  0.11157  1.31695    ignored, due to not correspond
        8  0.00308  0.00995  0.00669
        1 -1.16558 -0.13113  1.64921
        1  1.30697  1.80453  1.36739
    """

    fp = tempfile.NamedTemporaryFile()
    fp.write(txtfile.encode('utf-8'))
    # reset fp pointer
    fp.seek(0)
    rf = ReadFile(fp.name,debug=True)
    # fp will be destoried inside open!
    if not rf.nice:
        print('Error: on ReadFile')
        print(rf.info)
        return
    rf.run()

    sf = SaveFile(rf.system,fname='mytest.xyz',ftype='xsf')
    assert sf.nice
    assert sf.fname == 'mytest.xyz.xsf'
    assert sf.ftype == 'xsf'


    sf1 = SaveFile(rf.system,energy=rf.energy)
    if not sf1.nice:
        print('Error: on SaveFile')
        print(sf1.info)
        return
    sf1.run()


    mylist = [
        [
            [6, -0.02976,  0.02692,  1.34829],
            [6,  0.37245,  1.46435,  1.85092],
        ],
        [
            [8,  0.04793, -0.07239,  0.03349],
        ],
        [
            [7,  0.70929,  1.43837,  3.32682],
            [1, -0.40465,  2.19862,  1.65751],
            [6,  0.37245,  1.46435,  1.85092],
        ]
    ]
    sf2 = SaveFile(mylist)
    if not sf2.nice:
        print('Error: on SaveFile')
        print(sf2.info)
        return
    sf2.run()




class BondPerception:
    """Bond connecting perception

    Reference:
        Zhang, Q., et al.
        A rule-based algorithm for automatic bond type perception.
        J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26

    Args:
        system : single: 2D List[[atomtype, x,y,z], ...]
        userinputs (bool): if it is True, index in list starts at 1

    Attributes:
        conb        : all bonds connections
        bcon        : perceptive bonds connections
        nconb       : all nonbonds connections
        fconb       : fragments all bonds connections
        fnconb      : fragments nonbonds connections
        cfnconb     : cross fragments nonbonds connections
        fragments   : 2D List[ List[int], ... ]
        atradius    : atom radius   :   1D List[float]
    """
    def __init__(self,system,userinputs=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.system = system if isinstance(system,list) else []
        self.userinputs = True if userinputs is True else False

        if len(self.system) <= 2:
            self.nice = False
            self.info = 'Fatal: too less'
            return
        
        # get real atomtype list
        # format: 2D str: [ [sign, number, name], ... ]
        atominfo = [[i[0], str(i[1]), i[4]] for i in FAI.atominfo]
        self.atradius = []
        for atom in self.system:
            bo = True
            if not isinstance(atom,list) or len(atom) != 4:
                bo = False
            else:
                bo = False
                atype = str(atom[0]) if isinstance(atom[0],int) else atom[0].capitalize()
                for cnt,ndx in enumerate(atominfo):
                    if atype in ndx:
                        bo = True
                        r = FAI.atominfo[cnt][2]
                        break
            if not bo:
                self.nice = False
                self.info = 'Fatal: wrong defined: atom: {:}'.format(at[0])
                return
            self.atradius.append(r)



    def run(self):
        self.conb, self.bcon, self.nconb = self.calc_bcons(self.system,self.atradius)

        # based on self.bcon, calculate fragments
        # note: return is zero-based
        self.fragments = self.calc_bfs(self.bcon,len(self.system))

        self.fconb = []
        self.fnconb = []
        for ref in self.fragments:
            fc, fbc, fnbc = self.calc_bcons(self.system,self.atradius,reflist=ref)
            self.fconb.extend(fc)
            self.fnconb.extend(fnbc)
        
        self.cfnconb = []
        for i,ref in enumerate(self.fragments):
            j = i + 1
            while j < len(self.fragments):
                for t in ref:
                    for k in self.fragments[j]:
                        if t > k: t, k = k, t
                        if [t,k] not in self.cfnconb:
                            self.cfnconb.append([t,k])
                j += 1

        if self.userinputs:
            self.conb = [[i+1 for i in t] for t in self.conb]
            self.bcon = [[i+1 for i in t] for t in self.bcon]
            self.nconb = [[i+1 for i in t] for t in self.nconb]
            self.fragments = [[i+1 for i in t] for t in self.fragments]
            self.fconb = [[i+1 for i in t] for t in self.fconb]
            self.fnconb = [[i+1 for i in t] for t in self.fnconb]
            self.cfnconb = [[i+1 for i in t] for t in self.cfnconb]



    def calc_bcons(self,system,atradius,reflist=None):
        con = []
        bcon = []
        nconb = []
        reflist = list(range(len(system))) if reflist is None else reflist
        for i,ref in enumerate(system[:-1]):
            if i not in reflist: continue
            j = i + 1
            while j < len(system):
                if j in reflist:
                    atom = system[j]
                    dx = ref[1] - atom[1]
                    dy = ref[2] - atom[2]
                    dz = ref[3] - atom[3]
                    ds = dx*dx + dy*dy + dz*dz
                    rij = atradius[i] + atradius[j] + 0.4
                    if ds >= 0.64 and ds <= rij*rij:
                        bcon.append([i,j])
                    else:
                        nconb.append([i,j])
                    con.append([i,j])
                j += 1
        return con, bcon, nconb



    def calc_bfs(self,bcon,tot=None):
        """Breadth First Search to calculate 2D collections

        Return:
            bfs : 2D [fragments, fragments, ...]
        """
        if tot is None:
            tot = max([max(i) for i in bcon]) + 1
        visited = [False for i in range(tot)]
        bfs = []
        for ref in bcon:
            if visited[ref[0]]: continue
            queue = [ref[0]]
            visited[ref[0]] = True
            ls = []
            while queue:
                s = queue.pop()
                ls.append(s)
                for xmp in bcon:
                    if s == xmp[0] and not visited[xmp[1]]:
                        visited[xmp[1]] = True
                        queue.append(xmp[1])
                    if s == xmp[1] and not visited[xmp[0]]:
                        visited[xmp[0]] = True
                        queue.append(xmp[0])
            bfs.append(sorted(ls))
        flatten = []
        for i in bfs: flatten.extend(i)
        for i in range(tot):
            if i not in flatten:
                bfs.append([i])
        return bfs




class AnglePerception(BondPerception):
    """
    Attributes: (new)
        cona        : all angles connections
        acon        : perceptive angles connections
        nacon       : all nonangles connections
        fcona       : fragments all angles connections
        fncona      : fragments nonangles connections
        cfnacon     : cross two fragments nonangles connections
        cfncona     : cross three fragments nonangles connections
    """
    def __init__(self,system,*args,**kwargs):
        super().__init__(system,*args,**kwargs)
        if not self.nice: return



    def run(self):
        super().run()
        self.cona = []
        self.acon = []
        self.ncona = []
        self.fcona = []
        self.fncona = []
        self.cfnacon = []
        self.cfncona = []
        if len(self.system) < 3: return

        # take care results from BondPerception
        if self.userinputs:
            mybcon = [[i-1 for i in j] for j in self.bcon]
            myfragments = [[i-1 for i in j] for j in self.fragments]

        self.cona = self.calc_conas(list(range(len(self.system))))
        self.acon = self.calc_acons(mybcon)

        # deep copy
        for ref in self.cona:
            if ref not in self.acon:
                self.ncona.append([i for i in ref])

        for ref in myfragments:
            newfcona = self.calc_conas(ref)
            self.fcona.extend(newfcona)
            gfbcon = []
            for i in ref:
                for zmp in mybcon:
                    if i in zmp:
                        gfbcon.append(zmp)
            # it may have repeats
            newfbcon = []
            for i,ref in enumerate(gfbcon):
                j = i + 1
                bo = True
                while j < len(gfbcon):
                    if ref[0] in gfbcon[j] and ref[1] in gfbcon[j]:
                        bo = False
                        break
                    j += 1
                if bo:
                    newfbcon.append(ref)

            newfacon = self.calc_acons(newfbcon)
            # deep copy
            for ndx in newfcona:
                if ndx not in newfacon:
                    self.fncona.append([i for i in ndx])

        for i,ref in enumerate(myfragments):
            for j in range(i+1, len(myfragments)):
                for k in range(j+1, len(myfragments)):
                    for ai in ref:
                        for aj in myfragments[j]:
                            for ak in myfragments[k]:
                                # center is ai
                                if aj < ak:
                                    self.cfncona.append([aj,ai,ak])
                                else:
                                    self.cfncona.append([ak,ai,aj])
                                # center is aj
                                if ai < ak:
                                    self.cfncona.append([ai,aj,ak])
                                else:
                                    self.cfncona.append([ak,aj,ai])
                                # center is ak
                                if ai < aj:
                                    self.cfncona.append([ai,ak,aj])
                                else:
                                    self.cfncona.append([aj,ak,ai])

        for i,ref in enumerate(myfragments):
            if len(ref) < 2: continue
            for j in range(i+1, len(myfragments)):
                for s in range(len(ref)):
                    ai = ref[s]
                    for t in range(s+1,len(ref)):
                        aj = ref[t]
                        for ak in myfragments[j]:
                            # center is ai
                            if aj < ak:
                                self.cfnacon.append([aj,ai,ak])
                            else:
                                self.cfnacon.append([ak,ai,aj])
                            # center is aj
                            if ai < ak:
                                self.cfnacon.append([ai,aj,ak])
                            else:
                                self.cfnacon.append([ak,aj,ai])
                            # center is ak
                            if ai < aj:
                                self.cfnacon.append([ai,ak,aj])
                            else:
                                self.cfnacon.append([aj,ak,ai])

        if self.userinputs:
            self.cona = [[i+1 for i in j] for j in self.cona]
            self.acon = [[i+1 for i in j] for j in self.acon]
            self.ncona = [[i+1 for i in j] for j in self.ncona]
            self.fcona = [[i+1 for i in j] for j in self.fcona]
            self.fncona = [[i+1 for i in j] for j in self.fncona]
            self.cfnacon = [[i+1 for i in j] for j in self.cfnacon]
            self.cfncona = [[i+1 for i in j] for j in self.cfncona]



    def calc_conas(self,nlist):
        """calc all angles connections based on given list
        Args:
            nlist: 1D List[int]
        """
        if len(nlist) < 3: return []
        cona = []
        i = 0
        while i < len(nlist):
            j = i + 1
            while j < len(nlist):
                k = j + 1
                while k < len(nlist):
                    # center index i
                    if nlist[j] < nlist[k]:
                        cona.append([nlist[j], nlist[i], nlist[k]])
                    else:
                        cona.append([nlist[k], nlist[i], nlist[j]])
                    # center index j
                    if nlist[i] < nlist[k]:
                        cona.append([nlist[i], nlist[j], nlist[k]])
                    else:
                        cona.append([nlist[k], nlist[j], nlist[i]])
                    # center index k
                    if nlist[j] < nlist[i]:
                        cona.append([nlist[j], nlist[k], nlist[i]])
                    else:
                        cona.append([nlist[i], nlist[k], nlist[j]])
                    k += 1
                j += 1
            i += 1
        return cona



    def calc_acons(self,bcon):
        if len(bcon) < 2: return []
        # rule: acon: [i,j,k]   :   i < k
        acon = []
        for cnt,ref in enumerate(bcon):
            # center is ref[1]
            i, j = ref[0], ref[1]
            for num in range(len(bcon)):
                if num == cnt: continue
                if j == bcon[num][0]:
                    k = bcon[num][1]
                elif j == bcon[num][1]:
                    k = bcon[num][0]
                else:
                    k = None
                if k is not None:
                    if k < i: i, k = k, i
                    if [i,j,k] not in acon:
                        acon.append([i,j,k])
            # center is ref[0]
            j, i = ref[0], ref[1]
            for num in range(len(bcon)):
                if num == cnt: continue
                if j == bcon[num][0]:
                    k = bcon[num][1]
                elif j == bcon[num][1]:
                    k = bcon[num][0]
                else:
                    k = None
                if k is not None:
                    if k < i: i, k = k, i
                    if [i,j,k] not in acon:
                        acon.append([i,j,k])
        return acon




def test_class_Perception():
    def checkrepeats(mylist):
        for i,ref in enumerate(mylist):
            bo = False
            for j in range(i+1,len(mylist)):
                bo = True
                for k in range(len(ref)):
                    if ref[k] != mylist[j][k]:
                        bo = False
                        break
                if bo:
                    break
            if bo:
                return True
        return False


    def checklist(totlist,sublist,adjlist=None):
        """
        inputs: List[[int]]

        check: totlist == sublist + adjlist ?

        Return:
            True    :   equals
            False   :   wrong
            None    :   contains
        """
        if checkrepeats(totlist):
            print('Fatal: repeats on totlist')
            return False
        if checkrepeats(sublist):
            print('Fatal: repeats on sublist')
            return False

        if adjlist is None:
            adjlist = []
            if len(totlist) < len(sublist): return False
        else:
            if checkrepeats(adjlist):
                print('Fatal: repeats on adjlist')
                return False
            if len(totlist) < len(sublist) + len(adjlist): return False
        
        # check inner-length, do not modify the reference
        sublist = [i for i in sublist]
        sublist.extend(adjlist)

        tl = [len(i) for i in totlist]
        sl = [len(i) for i in sublist]
        if tl[0] != sl[0] or len(set(tl)) != 1 or len(set(sl)) != 1:
            return False

        chktlist = [None for i in range(len(totlist))]
        chkslist = [False for i in range(len(sublist))]
        for i,ref in enumerate(totlist):
            bo = False
            for j,xmp in enumerate(sublist):
                if chkslist[j]: continue
                bo = True
                for n in range(len(xmp)):
                    if ref[n] != xmp[n]:
                        bo = False
                        break
                if bo:
                    chkslist[j] = True
                    break
            if bo:
                chktlist[i] = True
        
        if False in chkslist:
            return False
        if None in chktlist:
            return None
        return True
    

    def checkobj(obj,info='mol'):
        errinfo = 'Error: {:}:'.format(info)
        print('Note: {:}: fragments='.format(info),obj.fragments)
        if checklist(obj.conb, obj.bcon, obj.nconb) != True:
            print(errinfo,'conb != bcon + nconb')
        if min([len(i) for i in obj.fragments]) >= 2:
            if checklist(obj.nconb, obj.fnconb, obj.cfnconb) != True:
                print(errinfo,'nconb != fnconb + cfnconb')
        else:
            if checklist(obj.nconb, obj.fnconb, obj.cfnconb) == False:
                print(errinfo,'nconb ?? fnconb + cfnconb')
        if checklist(obj.fconb, obj.bcon, obj.fnconb) != True:
            print(errinfo,'fconb != bcon + fnconb')
        if checklist(obj.cona, obj.acon, obj.ncona) == False:
            print(errinfo,'cona != acon + ncona')
        if len(obj.fragments) >= 3:
            if min([len(i) for i in obj.fragments]) <= 1:
                mycfa = [i for i in obj.fncona]
                mycfa.extend(obj.cfnacon)
                if checklist(obj.ncona, obj.fncona, mycfa) != True:
                    print(errinfo,'ncona != fncona + cfncona')
            else:
                if checklist(obj.ncona, obj.fncona, obj.cfncona) == False:
                    print(errinfo,'-1-: ncona ?? fncona + cfncona')
        else:
            if checklist(obj.ncona, obj.fncona, obj.cfncona) == False:
                    print(errinfo,'-2-: ncona ?? fncona + cfncona')
        if checklist(obj.fcona, obj.acon, obj.fncona) == False:
            print(errinfo,'fcona < acon + fncona')



    amol = [
        ['C',   -1.7051262,  -0.4420827,   1.9222438],
        ['H',   -1.3484718,  -1.4508927,   1.9222438],
        ['H',   -1.3484534,   0.0623155,   2.7958953],
        ['H',   -2.7751262,  -0.4420695,   1.9222438],
        ['C',   -1.1917840,   0.2838736,   0.6648389],
        ['H',   -1.5482973,  -0.2206384,  -0.2088121],
        ['H',   -1.5485980,   1.2926271,   0.6647411],
        ['H',   -0.1217840,   0.2840299,   0.6649364],
    ]
    #      4   6
    #      |   |
    #   3--1---5--7
    #      |   |
    #      2   8
    adict = {
        'userinputs': True,
        'fragments' : [[1,2,3,4,5,6,7,8]],
        'bcon'      : [[1,2],[1,3],[1,4],[1,5],[5,6],[5,7],[5,8]],
    }
    mf = AnglePerception(amol,**adict)
    mf.run()
    checkobj(mf,'amol')


    bmol = [
        ['C',   -1.7051262,  -0.4420827,   1.9222438],
        ['H',   -1.3484718,  -1.4508927,   1.9222438],
        ['H',   -1.3484534,   0.0623155,   2.7958953],
        ['H',   -2.7751262,  -0.4420695,   1.9222438],
        ['N',   -0.6491598,   2.6047218,  -1.2667619],
        ['H',   -0.3158379,   1.6619088,  -1.2667619],
        ['H',   -0.3158207,   3.0761220,  -0.4502651],
        ['H',   -0.3158207,   3.0761220,  -2.0832586],
    ]
    #      4           6
    #      |           |
    #   3--1--2     8--5--7
    bdict = {
        'userinputs': True,
        'fragments' : [[1,2,3,4],[5,6,7,8]],
    }
    mf = AnglePerception(bmol,**bdict)
    mf.run()
    checkobj(mf,'bmol')


    cmol = [
        ['C',   -1.7051262,  -0.4420827,   1.9222438],
        ['H',   -2.7751262,  -0.4420695,   1.9222438],
        ['N',   -0.6491598,   2.6047218,  -1.2667619],
        ['H',   -0.3158379,   1.6619088,  -1.2667619],
        ['H',   -0.3158207,   3.0761220,  -0.4502651],
        ['H',   -0.3158207,   3.0761220,  -2.0832586],
        ['O',    1.8267490,   4.3088580,  -2.2150034],
        ['H',    2.7867490,   4.3088580,  -2.2150034],
        ['H',    1.5062945,   5.2137938,  -2.2150034],
    ]
    #              5
    #              |
    #   1--2    4--3--6     8--7--9
    cdict = {
        'userinputs': True,
        'fragments' : [[1,2],[3,4,5,6],[7,8,9]]
    }
    mf = AnglePerception(cmol,**cdict)
    mf.run()
    checkobj(mf,'cmol')


    dmol = [
        ['C',  -2.8322784,   0.1898734,   0.0000000],
        ['H',  -2.4756239,  -0.8189366,   0.0000000],
        ['H',  -2.4756055,   0.6942716,   0.8736515],
        ['F',  -2.3822706,   0.8262637,  -1.1022706],
        ['N',   0.0203601,   1.2966940,   2.2276179],
        ['Cl',  0.5837031,   2.0933603,   3.6074974],
        ['H',   0.3536820,   0.3538809,   2.2276179],
        ['H',   0.3536992,   1.7680942,   1.4111212],
        ['O',   2.8367428,   2.6674334,   5.2572880],
        ['H',   2.5162882,   3.5723693,   5.2572880],
        ['H',   3.7967428,   2.6674334,   5.2572880],
    ]
    #      3          7
    #      |          |
    #   2--1--4    6--5--8     10--9--11
    ddict = {
        'userinputs': True,
        'fragments' : [[1,2,3,4],[5,6,7,8],[9,10,11]]
    }
    mf = AnglePerception(dmol,**ddict)
    mf.run()
    checkobj(mf,'dmol')




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
    def __init__(self,system,bcon=None,acon=None,btol=None,atol=None,
                obonds=None,oangles=None,opar=None,oall=None,
    
    *args,**kwargs):
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



    def get_reflist(self):
        """note: it is valid after self.run() is called"""
        return self.reflist



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

        self.reflist = self.get_reflist()

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
            # note: self.bondlist & self.anglelist are from parent
            for i in self.bondlist: bltot.append(i)
            for i in self.anglelist: altot.append(i)

            # 3rd, calculate reflist
            reflist = self.calc_mols_filter_list(bltot,altot,[binc,ainc])
            reflist = [i-ndxlth for i in reflist if i > ndxlth]
            self.reflist = reflist

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
    if fname is None: fname = 'filtration-image.png'

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
        
        if len(systemlist) == 0:
            print('Warning: no inputs')
            return

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

        # get fname & ftype
        FD = SaveFile([0],*self.args,**self.kwargs)
        self.kwargs['fname'] = file_gen_new(FD.fname,fextend=FD.ftype)

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
                fgp = file_gen_new(ftmp,fextend='png',foriginal=False)
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
            fga = file_gen_new(fga,fextend='png',foriginal=False)
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
            if len(DATA.system) != 0:
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
            CFF.run()

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
        help='Output file type, [txt, xsf, xyz]',
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
        'ftype'                 :   None,
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

