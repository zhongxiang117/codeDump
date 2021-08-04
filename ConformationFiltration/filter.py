#!/usr/bin/env python3

import os
import sys
import math
import argparse
import matplotlib.pyplot as plt
import tempfile
import random
import time


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
    'version 3.3.7  : refine BulkProcess',
    'version 3.3.8  : refine parsecmd',
    'version 3.4.0  : RELEASE',
    'version 3.5.0  : deep refine Filtration, RELEASE',
    'version 3.6.0  : add class PlotSamples',
    'version 3.7.0  : add subcommand plot, RELEASE',
    'version 3.7.1  : small fix on userinputs',
    'version 3.8.0  : add plot random seed',
    'version 3.9.0  : add dynamic and static method in Filtration',
    'version 4.0.0  : RELEASE',
    'version 4.0.1  : add more calculation type in Check',
    'version 4.1.0  : add report of execution time',
    'version 4.2.0  : prompt check on memory usage overhead',
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
        self.system = []
        self.energy = []
        self.debug = True if debug is True else False

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
            self.info = 'Fatal: file not support: {:}'.format(self.file)
            return

    def run(self):
        """process input file

        require:
            atomtype entry can be mixed, which means,

            a) 1  x  y  z
            b) H  x  y  z

            they are equivalent
        """
        if self.debug: print('Note: reading data from file: {:}'.format(self.file))
        prolist,enelist,errlist = getattr(self,'read_'+self.ext)()
        if not len(prolist): return

        if self.debug:
            for i in errlist:
                print('Warning: ignoring: {:}: {:}'.format(i[0],i[1]))

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
        for n,mol in enumerate(prolist):
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
                self.energy.append(enelist[n])
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
            sub = profile[i].strip()
            if not len(sub):
                i += 1
                continue
            if sub[0] == '#':
                ls = [[sub,i], ]
                j = i + 1
                while j < len(profile):
                    sub = profile[j].strip()
                    if not len(sub):
                        j += 1
                        continue
                    if sub[0] == '#': break
                    ls.append([sub,j])
                    j += 1
                promol.append(ls)
                i = j
            else:
                i += 1

        prolist = []
        enelist = []
        errlist = []
        for mol in promol:
            bo = False
            if len(mol) <= 2:
                bo = True
                errnum = mol[0][1] + 1
                errline = 'Wrong format'
            else:
                if mol[1][0] != 'ATOMS':
                    bo = True
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
                    else:
                        bo = True
                    if bo:
                        errnum = t[1] + 1
                        errline = t[0]
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
            if not len(sub):
                if not len(mol): continue
                promol.append(mol)
                # initialize
                mol = []
            else:
                mol.append([sub,cnt])
        # last mol
        if len(mol): promol.append(mol)

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
            if len(sub):
                mol.append([sub,cnt])
            else:
                if not len(mol): continue
                promol.append(mol)
                # initialize
                mol = []
        # last mol
        if len(mol): promol.append(mol)

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

        if not len(self.system):
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
            if self.fname[ndx:] != '.' + self.ftype:
                self.fname = self.fname[:ndx] + '.' + self.ftype

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
            if len(self.energy) and self.energy[ndx] is not None:
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
            if len(self.energy) and self.energy[ndx] is not None:
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
            if len(self.energy) and self.energy[ndx] is not None:
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
                self.info = 'Fatal: wrong defined: atom: {:}'.format(atom[0])
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
        mybcon = self.bcon
        myfragments = self.fragments
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
    """Filter molecules based on bonds & angles connections

    Inputs:
        system : 3D List[ List[[atomtype, x,y,z], ...], ...]
        userinputs (bool): if it is True, index in list starts at 1

        bcon : System[[atom-i, atom-j], ...] : List[[int,int], ...]
        acon : System[[atom-i, atom-j, atom-k], ...] : List[[int,int,int], ...]

        btol : tolerance on bond, Angstrom
        atol : tolerance on angle, degree

        obpar  :  Boolean  :  whether calculate bonds prob_bpar  :  default False
        oball  :  Boolean  :  whether calculate bonds prob_ball  :  default True
        oapar  :  Boolean  :  whether calculate angles prob_apar  :  default False
        oaall  :  Boolean  :  whether calculate angles prob_aall  :  default True


    Attributes:
        system  :  good molecules after filtration
        sysbad  :  filtered out molecules

        prob_begin : begin probability  :   dict
            keys:
                bpar    : 3D : List[ [List[int],float] ]
                ball    : 2D : List[ List[int],  float ]
                apar    : 3D : List[ [List[int],float] ]
                aall    : 2D : List[ List[int],  float ]

        prob_final :  same format as prob_begin

        bondlist   :  2D  :  List[ List[float] ]  :  good, correspond to bcon
        anglelist  :  2D  :  List[ List[float] ]  :  good, correspond to acon

        reflist    :  1D  List[int]     : sorted index of for bad molecules


    Caution:
            Because the BOSS solute move is only for the single time move,
            so Sum(difference for untouched atoms) =~ 0.0,
            thus, when doing filtration, please take care of btol & atol.
    """
    def __init__(self,system=None,keepndxlist=None,userinputs=None,
                bcon=None,acon=None,btol=None,atol=None,seed=None,
                mode=None,vndx=None,borandom=None,boall=None,
                obpar=None,oball=None,oapar=None,oaall=None,
                *args,**kwargs):
        self.system = system
        self.keepndxlist = keepndxlist
        self.userinputs = True if userinputs is True else False

        self.bcon = [] if bcon is None else bcon
        self.acon = [] if acon is None else acon
        self.btol = 0.1 if btol is None else btol   # Angstrom
        self.atol = 0.1 if atol is None else atol   # degree

        self.seed = seed if seed else random.randrange(100000000)
        random.seed(self.seed)

        if mode is None or mode.lower() in ['d','dynamic']:
            self.mode = 'dynamic'
        else:
            self.mode = 'static'
        self.vndx = vndx
        self.borandom = borandom
        self.boall = True if boall is None else boall

        self.obpar = True if obpar is True else False
        self.oball = False if oball is False else True
        self.oapar = True if oapar is True else False
        self.oaall = False if oaall is False else True

        if self.userinputs:
            self.userinputs = False
            self.bcon = [[i-1 for i in j] for j in self.bcon]
            self.acon = [[i-1 for i in j] for j in self.acon]
        self.prob_begin = {'bpar':[], 'ball':[], 'apar':[], 'aall':[]}
        self.prob_final = {'bpar':[], 'ball':[], 'apar':[], 'aall':[]}

    def run(self):
        """attemption on filtering
        """
        # to improve efficiency, bondlist only needs to be calculated once
        if len(self.bcon): print('Note: calculating bonds connections ...')
        bondlist = self.calc_square_distance(self.system,self.bcon)
        if len(self.acon): print('Note: calculating angles connections ...')
        anglelist = self.calc_angle_degree(self.system,self.acon)

        # increments
        binc = self.btol * self.btol
        ainc = self.atol

        if self.keepndxlist is None: self.keepndxlist = []
        # important, increase efficiency
        self.keepndxlist = sorted(self.keepndxlist)

        if self.obpar or self.oball:
            print('Note: calculating begin bonds probability ...')
            self.prob_begin['bpar'], self.prob_begin['ball'] = self.calc_probs(bondlist,binc,self.obpar,self.oball)
        if self.oapar or self.oaall:
            print('Note: calculating begin angles probability ...')
            self.prob_begin['apar'], self.prob_begin['aall'] = self.calc_probs(anglelist,ainc,self.oapar,self.oaall)
        print('Note: calculating repeats reference ...')
        self.reflist = self.calc_filterlists(bondlist,anglelist,binc,ainc,mode=self.mode,vndx=self.vndx,
                                            borandom=self.borandom,boall=self.boall,keepndxlist=self.keepndxlist)
        
        print('Note: updating ...')
        tmpsys = []
        self.bondlist = []
        self.anglelist = []
        self.sysbad = []
        cnt = 0
        trn = 0
        self.reflist.append(-1)
        self.keepndxlist.append(-1)
        for ndx in range(len(self.system)):
            if ndx == self.reflist[cnt]:
                cnt += 1
                self.sysbad.append(self.system[ndx])
            else:
                if ndx == self.keepndxlist[trn]:
                    trn += 1
                else:
                    tmpsys.append(self.system[ndx])
                    if len(bondlist): self.bondlist.append(bondlist[ndx])
                    if len(anglelist): self.anglelist.append(anglelist[ndx])
        self.fratio = 1.0 - len(tmpsys)/len(self.system)
        # alias
        self.system = tmpsys
        self.reflist.pop(len(self.reflist)-1)
        self.keepndxlist.pop(len(self.keepndxlist)-1)
        
        if self.obpar or self.oball:
            print('Note: calculating final bonds probability ...')
            self.prob_final['bpar'], self.prob_final['ball'] = self.calc_probs(self.bondlist,binc,self.obpar,self.oball)
        if self.oapar or self.oaall:
            print('Note: calculating final angles probability ...')
            self.prob_final['apar'], self.prob_final['aall'] = self.calc_probs(self.anglelist,ainc,self.oapar,self.oaall)

    def calc_probs(self,datalist,dt,opar=None,oall=None):
        """
        Inputs:
            dt : float

        Return:
            prob_par  : 3D : List[ [List[int],float] ]

                => the probability of any two molecules whose difference fall
                   within defined range

                explanation:
                    assume three molecules, the number of atoms are 4;
                    A = [a1, a2, a3, a4]
                    B = [b1, b2, b3, b4]
                    C = [c1, c2, c3, b4]

                    define its connection index is: [[0,1], [1,2], [1,3]]
                    means bonded atom pairs in A are [a1-a2, a2-a3, a2-a4]

                    now, calculate their bondlist;

                    DA = [Da01, Da12, Da13]
                    DB = [Db01, Db12, Db13]
                    DC = [Dc01, Dc12, Dc13]

                    thus, the difference of bond-length for any two molecules
                    corresponding to connection will be DA-DB, DA-DC, DB-DC

                    therefore a new array will be got:
                    New = [Da01-Db01, Da01-Dc01, Db01-Dc01]

                    finally we can calculate molecule-increments probability


            prob_all  : 2D : List[ List[int],  float ]

                => means probability on overall difference

                explanation:
                    continuing top, we have already calculated DA, DB, DC

                    calculate their average values;

                    TA = Average(DA)
                    TB = Average(DB)
                    TC = Average(DC)

                    following the same rule, calculate differences,
                    New = [TA-TB, TA-TC, TB-TC]

                    finally we can calculate overall probability

        Note:
            results are not normalized
        """
        def calc_list(ls,dt):
            stls = sorted(ls)
            rmin = stls[0]
            steps = (stls[-1] - rmin) / dt
            # for round-of-errors, as well as endpoint
            steps = int(steps) + 2
            # all on same values
            if steps <= 3:
                return [len(ls)],rmin
            prolist = []
            ndx = 0
            # for endpoint
            stls.append(stls[-1]+dt)
            for v in range(1,steps):
                high = rmin + v*dt
                for n,t in enumerate(stls[ndx:]):
                    if t >= high:
                        prolist.append(n)
                        ndx += n
                        break
            return prolist,rmin

        if len(datalist) <= 1: return [],[]

        prob_par = []
        if opar:
            print('    --> computing on par entry ...')
            for cnt in range(len(datalist[0])):
                ls = [t[cnt] for t in datalist]
                p = calc_list(ls,dt)
                prob_par.append(p)
        prob_all = []
        if oall:
            print('    --> computing on system ...')
            # Caution! only one movement
            #sub = [sum(i)/len(i) for i in datalist]
            sub = [sum(i) for i in datalist]
            prob_all = calc_list(sub,dt*len(datalist[0]))
        return prob_par, prob_all

    def calc_filterlists(self,bondlist,anglelist,binc,ainc,mode=None,vndx=None,
                        borandom=None,boall=None,keepndxlist=None):
        """
        Rule:
            Since the variance is only based on the single movement,
            either can be a bond or angle, so the changes on any two adjacent
            molecules are smaller than tolerance will be removed,
            thus, dt should not multiple the total numbers,
            because Sum(mol untouched atoms) =~ 0.0

        Return:
            reflist  :  List[int]  :  index of molecules waiting to be removed
        """
        if len(bondlist) <= 3 and len(anglelist) <= 3: return []
        bl = [sum(i) for i in bondlist]
        al = [sum(i) for i in anglelist]
        if mode is None or mode.lower() in ['dynamic','d']:
            if boall:
                return self._calc_filterlists_all_dynamic(bl,al,binc,ainc,keepndxlist)
            # dynamic-separate on bondlist
            reflist = self._calc_filterlists_sep_dynamic(bl,al,binc,ainc,keepndxlist)
            # dynamic-separate on anglelist
            reflist.append(-1)
            cnt = 0
            mybl = []
            myal = []
            ndxlist = []
            for i in range(len(bl)):
                if i == reflist[cnt]:
                    cnt += 1
                else:
                    mybl.append(bl[i])
                    myal.append(al[i])
                    ndxlist.append(i)
            reflist.pop(len(reflist)-1)
            tmplist = self._calc_filterlists_sep_dynamic(myal,mybl,ainc,binc,keepndxlist)
            reflist.extend([ndxlist[i] for i in tmplist])
            return sorted(reflist)
        return self._calc_filterlists_static(bl,al,binc,ainc,keepndxlist,vndx,borandom)

    def _calc_filterlists_all_dynamic(self,bl,al,binc,ainc,keepndxlist):
        """This part does is to recursively remove index, whose difference
           separately with two adjacent values is smaller than tolerance
        
        for example, we have inputs like;
        
          bal      = [0.0, 0.1, 0.15, 0.18, 0.19, 0.20, 0.3, 0.4, 0.43, 0.5]
        ndxlist    =   0    1    2     3     4     5     6    7    8     9
        
        inc = 0.1
        
        then we will know if we remove index in [2,3,4,8], the new array
        balnew = [0.0, 0.1, 0.20, 0.3, 0.4, 0.5] will meet the requirement

        sort by index
        Caution: for future debug, bal is not sorted"""
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]
        keepndxlist.append(-1)
        
        if not len(bl) or not len(al):
            if not len(bl):
                bal = al
                inc = ainc
            else:
                bal = bl
                inc = binc
        else:
            bal = [v+al[i] for i,v in enumerate(bl)]
            inc = binc + ainc
        nlist = sorted(range(len(bal)),key=lambda k: bal[k])
        mdel = [bal[j] - bal[nlist[i]] for i,j in enumerate(nlist[1:])]

        reflist = []
        ndx = 0
        cnt = 0
        while ndx < len(mdel):
            if mdel[ndx] < inc:
                v1 = nlist.pop(ndx+1)
                if v1 == keepndxlist[cnt]:
                    cnt += 1
                    if ndx >= len(mdel)-1: break
                    # be aware in here, nlist has been processed
                    v2 = nlist[ndx+1]
                    mdel.pop(ndx)
                    mdel[ndx] = bal[v2] - bal[v1]
                else:
                    reflist.append(v1)
                    if ndx >= len(mdel)-1: break
                    dt = mdel.pop(ndx+1)
                    mdel[ndx] += dt
            else:
                ndx += 1
        # reflist has to be sorted from smaller to bigger for following refinement
        return sorted(reflist)
    
    def _calc_filterlists_sep_dynamic(self,bl,al,binc,ainc,keepndxlist):
        if not len(bl) or not len(al):
            return self._calc_filterlists_all_dynamic(bl,al,binc,ainc,keepndxlist)
        
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]
        keepndxlist.append(-1)

        nlist = sorted(range(len(bl)),key=lambda k: bl[k])
        mdel = [bl[j] - bl[nlist[i]] for i,j in enumerate(nlist[1:])]
        reflist = []
        ndx = 0
        cnt = 0
        while ndx < len(mdel):
            if mdel[ndx] < binc:
                v0 = nlist[ndx]
                v1 = nlist[ndx+1]
                # be aware of negative value
                if abs(al[v1]-al[v0]) < ainc:
                    nlist.pop(ndx+1)
                    if v1 == keepndxlist[cnt]:
                        cnt += 1
                        if ndx >= len(mdel)-1: break
                        v2 = nlist[ndx+1]
                        mdel.pop(ndx)
                        mdel[ndx] = bl[v2] - bl[v1]
                    else:
                        reflist.append(v1)
                        if ndx >= len(mdel)-1: break
                        dt = mdel.pop(ndx+1)
                        mdel[ndx] += dt
                else:
                    ndx += 1
            else:
                ndx += 1
        return sorted(reflist)

    def _calc_filterlists_static(self,bl,al,binc,ainc,keepndxlist,vndx=None,borandom=None):
        # make a deep copy of keepndxlist
        if keepndxlist is None: keepndxlist = []
        keepndxlist = [i for i in keepndxlist]

        if not len(bl) or not len(al):
            if not len(bl):
                bal = al
                inc = ainc
            else:
                bal = bl
                inc = binc
        else:
            bal = [v+al[i] for i,v in enumerate(bl)]
            inc = binc + ainc

        nlist = sorted(range(len(bal)),key=lambda k: bal[k])
        vlist = [bal[i] for i in nlist]

        # always make vndx one-inc less than smallest value
        if vndx is None:
            vndx = vlist[0]
        elif vndx > vlist[0]:
            while vndx > vlist[0]:
                vndx -= inc
        else:
            while vndx < vlist[0]:
                vndx += inc
            vndx -= inc
        
        reflist = []
        cnt = 0
        n = int((vlist[-1]-vndx)/inc) + 1
        tot = len(vlist)
        for i in range(1,n):
            t = vndx + i*inc
            ls = []
            while cnt < tot:
                if vlist[cnt] >= t: break
                ls.append(nlist[cnt])
                cnt += 1
            if len(ls):
                bo = False
                # important, increase efficiency
                if len(keepndxlist):
                    for k in ls:
                        if k in keepndxlist:
                            bo = True
                            keepndxlist.remove(k)
                if bo:
                    reflist.extend(ls)
                elif len(ls) >= 2:
                    if borandom:
                        x = random.randrange(len(ls))
                        ls.pop(x)
                        reflist.extend(ls)
                    else:
                        reflist.extend(ls[1:])
        return sorted(reflist)

    def calc_square_distance(self,system,bcon):
        """
        Return:
            bondlist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]
        """
        if not len(bcon): return []
        bondlist = []
        for mol in system:
            ls = []
            for ndx in bcon:
                at1 = mol[ndx[0]]
                at2 = mol[ndx[1]]
                dx = at1[1] - at2[1]
                dy = at1[2] - at2[2]
                dz = at1[3] - at2[3]
                tmp = dx*dx + dy*dy + dz*dz
                ls.append(tmp)
            bondlist.append(ls)
        return bondlist

    def calc_angle_degree(self,system,acon):
        """
        Rule:
            assume cooridnates,

            A (ax, ay, az)
            B (bx, by, bz)
            C (cx, cy, cz)

            Vector,

            BA = (ax-bx, ay-by, az-bz)
            BC = (cx-bx, cy-by, cz-bz)

            cos<ABC> = Sum(bai*bci) / dis(BA) * dis(BC)

            <ABC> = math.acos(sigma) * 180.0 / math.pi

        Return:
            anglelist : 2D : System[Mol[l1,l2, ...], ...] : List[List[float]]
        """
        if not len(acon): return []
        anglelist = []
        cvt = 180.0 / math.pi
        for mol in system:
            ls = []
            for ndx in acon:
                a = mol[ndx[0]]
                b = mol[ndx[1]]
                c = mol[ndx[2]]
                ba = [a[1]-b[1],a[2]-b[2],a[3]-b[3]]
                bc = [c[1]-b[1],c[2]-b[2],c[3]-b[3]]
                tot = ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]
                sub = sum([i*i for i in ba]) * sum([i*i for i in bc])
                rst = math.acos(tot/pow(sub,0.5)) * cvt
                ls.append(rst)
            anglelist.append(ls)
        return anglelist


def test_class_Filtration_dynamic():
    """
    Be aware of the testing data file is used
    """
    def checkequal(alist,blist,tol=None):
        if tol is None: tol = 10**(-5)
        if len(alist) != len(blist): return False
        for i,la in enumerate(alist):
            for t,va in enumerate(la):
                if va - blist[i][t] > tol:
                    return False
        return True

    rf = ReadFile('choosetest.txt')
    if not rf.nice:
        print(rf.info)
        exit()
    rf.run()
    ltmp = list(range(len(rf.system[0])))
    if len(ltmp) <= 4:
        print('Fatal: atom too less, cannot debug')
        exit()
    bcon = []
    for i in range(5):
        s = random.sample(ltmp,k=2)
        st = sorted(s)
        if st not in bcon: bcon.append(st)
    acon = []
    for i in range(5):
        s = random.sample(ltmp,k=3)
        st = sorted(s)
        if st not in acon: acon.append(st)
    
    # dynamic-separate
    fd = {
        'system'    :   rf.system,
        'userinputs':   False,
        'bcon'      :   bcon,
        'acon'      :   acon,
    }
    fn = Filtration(**fd)
    fn.run()

    rawbondlist = fn.calc_square_distance(rf.system,fn.bcon)
    rawanglelist = fn.calc_angle_degree(rf.system,fn.acon)
    probondlist = [v for i,v in enumerate(rawbondlist) if i not in fn.reflist]
    proanglelist = [v for i,v in enumerate(rawanglelist) if i not in fn.reflist]
    assert checkequal(probondlist, fn.bondlist)
    assert checkequal(proanglelist, fn.anglelist)

    bl = [sum(i) for i in probondlist]
    al = [sum(i) for i in proanglelist]
    nl = sorted(range(len(bl)),key=lambda k: bl[k])
    stbl = [bl[i] for i in nl]
    dtbl = [v-stbl[i] for i,v in enumerate(stbl[1:])]
    bocheckfn = []
    binc = fn.btol * fn.btol
    for i,t in enumerate(dtbl):
        bo = False
        if t >= binc:
            bo = True
        else:
            if abs(al[nl[i+1]]-al[nl[i]]) >= fn.atol: bo = True
        bocheckfn.append(bo)
    assert all(bocheckfn)
    print('==> dynamic-sep filtration ratio: ',fn.fratio)

    # test on fixed cross filtration
    n = len(rf.system)
    fd['keepndxlist'] = random.sample(range(n),random.randint(2,n))
    fv = Filtration(**fd)
    fv.run()
    bl = [sum(i) for i in fv.bondlist]
    al = [sum(i) for i in fv.anglelist]
    nl = sorted(range(len(bl)),key=lambda k: bl[k])
    stbl = [bl[i] for i in nl]
    dtbl = [v-stbl[i] for i,v in enumerate(stbl[1:])]
    binc = fv.btol * fv.btol
    for i,t in enumerate(dtbl):
        bo = False
        if t >= binc:
            bo = True
        else:
            if abs(al[nl[i+1]]-al[nl[i]]) >= fv.atol: bo = True
        assert bo

    # dynamic-all
    fd = {
        'system'    :   rf.system,
        'userinputs':   False,
        'bcon'      :   bcon,
        'acon'      :   acon,
        'boall'     :   True,
    }
    fn = Filtration(**fd)
    fn.run()

    rawbondlist = fn.calc_square_distance(rf.system,fn.bcon)
    rawanglelist = fn.calc_angle_degree(rf.system,fn.acon)
    probondlist = [v for i,v in enumerate(rawbondlist) if i not in fn.reflist]
    proanglelist = [v for i,v in enumerate(rawanglelist) if i not in fn.reflist]
    assert checkequal(probondlist, fn.bondlist)
    assert checkequal(proanglelist, fn.anglelist)

    rawtol = fn.atol + fn.btol*fn.btol
    prototal = [sum(v)+sum(proanglelist[i]) for i,v in enumerate(probondlist)]
    stprototal = sorted(prototal)
    prodel = [v-stprototal[i] for i,v in enumerate(stprototal[1:])]
    bocheckfn = [True if i >= rawtol else False for i in prodel]
    assert all(bocheckfn)
    print('==> dynamic-all filtration ratio: ',fn.fratio)

    # test on fixed cross filtration
    n = len(rf.system)
    #fd['keepndxlist'] = random.sample(range(n),random.randint(2,n))
    fv = Filtration(**fd)
    fv.run()
    profv = [sum(v)+sum(fv.anglelist[i]) for i,v in enumerate(fv.bondlist)]
    stprofv = sorted(profv)
    prodt = [v-stprofv[i] for i,v in enumerate(stprofv[1:])]
    bocheckfv = [True if i >= rawtol else False for i in prodt]
    assert all(bocheckfv)
    ratio = 1.0 - len(set([*fv.reflist,*fv.keepndxlist]))/len(rf.system)
    print('==> dynamic-all filtration ratio: --keepndxlist ',ratio)


def test_class_Filtration_static():
    """
    Be aware of the testing data file is used
    """
    def checkequal(alist,blist,tol=None):
        if tol is None: tol = 10**(-5)
        if len(alist) != len(blist): return False
        for i,la in enumerate(alist):
            for t,va in enumerate(la):
                if va - blist[i][t] > tol:
                    return False
        return True

    rf = ReadFile('choosetest.txt')
    if not rf.nice:
        print(rf.info)
        exit()
    rf.run()
    ltmp = list(range(len(rf.system[0])))
    if len(ltmp) <= 4:
        print('Fatal: atom too less, cannot debug')
        exit()
    bcon = []
    for i in range(5):
        s = random.sample(ltmp,k=2)
        st = sorted(s)
        if st not in bcon: bcon.append(st)
    acon = []
    for i in range(5):
        s = random.sample(ltmp,k=3)
        st = sorted(s)
        if st not in acon: acon.append(st)
    
    # static
    fd = {
        'system'    :   rf.system,
        'userinputs':   False,
        'bcon'      :   bcon,
        'acon'      :   acon,
        'boall'     :   True,
        'mode'      :   'static',
        'vndx'      :   None,
        'borandom'  :   None,
    }
    fn = Filtration(**fd)
    fn.run()

    rawbondlist = fn.calc_square_distance(rf.system,fn.bcon)
    rawanglelist = fn.calc_angle_degree(rf.system,fn.acon)
    probondlist = [v for i,v in enumerate(rawbondlist) if i not in fn.reflist]
    proanglelist = [v for i,v in enumerate(rawanglelist) if i not in fn.reflist]
    assert checkequal(probondlist, fn.bondlist)
    assert checkequal(proanglelist, fn.anglelist)
    print(fn.fratio)

    if len(fn.acon):
        tol = fn.btol * fn.btol
    elif len(fn.bcon):
        tol = fn.atol
    else:
        tol = fn.atol + fn.btol*fn.btol
    raw = [sum(v)+sum(rawanglelist[i]) for i,v in enumerate(rawbondlist)]
    rmin = min(raw)
    bal = [sum(v)+sum(proanglelist[i]) for i,v in enumerate(probondlist)]
    stbal = sorted(bal)
    bocheckfn = []
    cnt = 0
    for i in range(1,len(stbal)-1):
        high = rmin + i*tol
        if stbal[cnt] < high:
            cnt += 1
            if stbal[cnt] < high:
                assert False

    # test on fixed cross filtration
    n = len(rf.system)
    fd['keepndxlist'] = random.sample(range(n),random.randint(2,n))
    fv = Filtration(**fd)
    fv.run()
    fv.reflist.extend(fd['keepndxlist'])
    probondlist = [v for i,v in enumerate(rawbondlist) if i not in fv.reflist]
    proanglelist = [v for i,v in enumerate(rawanglelist) if i not in fv.reflist]
    raw = [sum(v)+sum(rawanglelist[i]) for i,v in enumerate(rawbondlist)]
    rmin = min(raw)
    bal = [sum(v)+sum(proanglelist[i]) for i,v in enumerate(probondlist)]
    stbal = sorted(bal)
    bocheckfn = []
    cnt = 0
    for i in range(1,len(stbal)-1):
        high = rmin + i*tol
        if stbal[cnt] < high:
            cnt += 1
            if stbal[cnt] < high:
                assert False


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


def plot_save_image(ini,fin=None,dt=None,fname=None,key=None):
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
    if len(ini[0]) <= 3: boi = False
    if len(fin[0]) <= 3: bof = False
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


def getrealsizeof(o):
    """recursively get the real size of built-in objects, unit in bytes
    """
    tot = sys.getsizeof(o)
    if isinstance(o,(int,str,float)):
        return tot
    if isinstance(o,(list,tuple)):
        return tot + sum([getrealsizeof(i) for i in o])
    if isinstance(o,dict):
        return tot + sum(getrealsizeof(k)+getrealsizeof(v) for k,v in o.items())
    return tot


class BulkProcess:
    """bulk process for datafilelist based on indexfilelist

    Inputs:
        bcon (str|list): conb, bcon, ...
        acon (str|list): cona, acon, ...

    Attributes:
        overall_system  :   final good system
        overall_sysbad  :   all bad system

        overall_prob_begin  :   based on read data file
        overall_prob_final  :   based on read data file

    Note:
        fragments is input as human-readable number, starting at 1
    """
    def __init__(self,datafilelist=None,indexfilelist=None,
                bool_force_double_check=None,*args,**kwargs):
        self.nice = True
        self.info = ''
        self.mytime = time.time()
        self.datafilelist = []
        if datafilelist is not None:
            for f in datafilelist:
                if os.path.isfile(f):
                    self.datafilelist.append(f)
                else:
                    print('Warning: not a data file < {:} >, ignoring'.format(f))

        if not len(self.datafilelist):
            self.nice = False
            self.info = 'Fatal: no inputs'
            return

        self.indexfilelist = []
        if indexfilelist is not None:
            for f in indexfilelist:
                if os.path.isfile(f):
                    self.indexfilelist.append(f)
                else:
                    print('Warning: not an index file < {:} >, ignoring'.format(f))

        self.bool_force_double_check = False if bool_force_double_check is False else True
        self.args = args
        self.kwargs = kwargs

    def run(self,debug=None):
        fsize = sum([os.stat(i).st_size for i in self.datafilelist])
        fsize += sum([os.stat(i).st_size for i in self.indexfilelist])
        # set warning of maximum valid file size
        memmax = 500
        if fsize/1024/1024 > memmax:
            rsize = 0
            msize = 0
            fsize = 0
            total = 100000
            cnt = 0
            for file in [*self.datafilelist, *self.indexfilelist]:
                if not (file.endswith('.xyz') or file.endswith('.xsf') or file.endswith('.txt')):
                    continue
                fsize += os.stat(file).st_size
                if cnt > total: continue
                profile = []
                with open(file,mode='rt') as f:
                    while True:
                        line = f.readline()
                        if not len(line): break
                        profile.append(line)
                        if cnt > total: break
                        cnt += 1
                ftmp = tempfile.NamedTemporaryFile()
                ftmp.write(''.join(profile).encode('utf-8'))
                ftmp.flush()
                rf = ReadFile(ftmp.name,ext=file[file.rfind('.')+1:],debug=False)
                if len(profile) < 3000: rf.nice = False
                if rf.nice:
                    rf.run()
                    if len(rf.system):
                        rsize += os.stat(ftmp.name).st_size
                        msize += getrealsizeof(rf.system) + getrealsizeof(rf.energy)
                ftmp.close()
            
            if fsize/1024/1024 > memmax and rsize > 1.0:
                # convert to MB, get the 1.2 times memory
                fsize = fsize / rsize * msize / 1000 / 1000 * 1.2
                print('Warning: your input is super large')
                print('Warning: memory requested roughly will be: {:} MB'.format(round(fsize,2)))
                print('Be sure that you want to continue? y/yes, else not. Input: ',end='')
                if input().lower() not in ['y','yes']:
                    print('Note: you decided to quit, nothing will be processed')
                    return
                print()

        systemlist,energylist = self.get_datalist(self.datafilelist)
        if not sum([len(i) for i in systemlist]):
            self.nice = False
            self.info = 'Fatal: no inputs after process'
            return
        self.molnms = [len(i) for i in systemlist]

        sysndxlist = []
        if len(self.indexfilelist):
            sysndxlist,tmp = self.get_datalist(self.indexfilelist)

        # connections only need to be calculated once
        choose = [i for i in systemlist if len(i)]
        self.get_connections(choose[0][0])
        if not self.nice: return

        # to make cross filtration happen, sysndxlist should be at the first
        allsystem = []
        allenergy = []
        for i in sysndxlist:
            allenergy.extend([None for j in range(len(i))])
            allsystem.extend(i)
        allkeeps = list(range(len(allsystem)))
        for i in systemlist: allsystem.extend(i)
        for i in energylist: allenergy.extend(i)
        mf = Filtration(system=allsystem,keepndxlist=allkeeps,*self.args,**self.kwargs)

        # prompt for double check
        if self.bool_force_double_check:
            print('\nCheck: current work path:')
            print('   => {:}'.format(os.path.abspath('.')))
            print('\nCheck: data files:')
            for cnt,fd in enumerate(self.datafilelist):
                print('   => {:} -- molnms {:}'.format(fd,self.molnms[cnt]))

            if len(self.indexfilelist):
                print('Check: index files:')
                for cnt,fd in enumerate(self.indexfilelist):
                    print('   => {:} -- molnms {:}'.format(fd,len(sysndxlist[cnt])))

            print('Check: molecule fragments:')
            lt = []
            for i in self.kwargs['fragments']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: bond connection:')
            lt = []
            for i in self.kwargs['bcon']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: angle connection:')
            lt = []
            for i in self.kwargs['acon']: lt.append([j+1 for j in i])
            print('   => {:}'.format(lt))

            print('Check: total inputs < {:} >'.format(sum(self.molnms)))

            stmp = 'ON' if mf.oball else 'OFF'
            print('Check: (image) bonds all probability < {:} >'.format(stmp))
            stmp = 'ON' if mf.obpar else 'OFF'
            print('Check: (images) bonds par probability < {:} > (time consuming)'.format(stmp))
            stmp = 'ON' if mf.oaall else 'OFF'
            print('Check: (image) angles all probability < {:} >'.format(stmp))
            stmp = 'ON' if mf.oapar else 'OFF'
            print('Check: (images) angles par probability < {:} > (time consuming)'.format(stmp))
            print('Check: bonds tolerance < {:} Angstrom >'.format(mf.btol))
            print('Check: angles tolerance < {:} degree >'.format(mf.atol))
            if mf.mode == 'dynamic':
                if mf.boall:
                    print('Check: calculation type: < dynamic/all >')
                else:
                    print('Check: calculation type: < dynamic/separate >')
            else:
                if mf.borandom:
                    print('Check: calculation type: < static/random >')
                else:
                    print('Check: calculation type: < static/lowest-bit >')
                if mf.vndx: print('Check: calculation vndx: < {:} >'.format(mf.vndx))
            imtot = 0
            if mf.oball: imtot += 1
            if mf.oaall: imtot += 1
            if mf.obpar: imtot += len(mf.bcon)
            if mf.oapar: imtot += len(mf.acon)
            print('Check: number of images will be generated: < {:} >'.format(imtot))

            print('\nDo you want to continue? y/yes, else not. Input: ',end='')
            if input().lower() not in ['y','yes']:
                print('Note: you decided to quit, nothing will be processed')
                return
            print()

        if debug: return allsystem
        mf.run()
        self.seed = mf.seed
        self.mode = mf.mode
        self.boall = mf.boall
        self.vndx = mf.vndx
        self.borandom = mf.borandom

        # accumulation only on datafilelist
        tot = 0
        acclist = []
        for i in sysndxlist: tot += len(i)
        acclist.append(tot)
        for i in systemlist:
            tot += len(i)
            acclist.append(tot)

        # test
        #mf.reflist = [i for i in acclist[:-1]]
        # expect: rmlist = [1 for i in range(len(acclist)-1)]
        # mf.reflist is the sorted list
        self.rmnmlist = [0 for i in range(len(acclist)-1)]
        ndx = 1
        i = 0
        while ndx < len(acclist):
            cnt = 0
            while i < len(mf.reflist):
                if mf.reflist[i] < acclist[ndx]:
                    cnt += 1
                    i += 1
                else:
                    break
            if cnt == 0:
                ndx += 1
            else:
                self.rmnmlist[ndx-1] = cnt

        self.overall_energy = []
        self.overall_system = []
        totreflist = sorted([*allkeeps,*mf.reflist])
        totreflist.append(-1)
        ndx = 0
        for i,v in enumerate(allenergy):
            if i == totreflist[ndx]:
                ndx += 1
            else:
                self.overall_energy.append(v)
                self.overall_system.append(allsystem[i])
        self.bcon = mf.bcon
        self.acon = mf.acon
        self.btol = mf.btol
        self.atol = mf.atol
        self.boim = True if mf.oball or mf.obpar or mf.oaall or mf.oapar else False
        self.overall_prob_begin = mf.prob_begin
        self.overall_prob_final = mf.prob_final
        self.save_files()

    def save_files(self):
        print('\nNote: saving bulk process results ...')
        tot = len(self.overall_system)
        print('Note: total molnms: < {:} >'.format(sum(self.molnms)))
        print('Note: final molnms: < {:} >'.format(tot))
        ratio = ('%f' % (1-round(tot/sum(self.molnms),6))).rstrip('0').rstrip('.')
        print('Note: filtration ratio: < {:} >'.format(ratio))

        # files
        fd = SaveFile(self.overall_system,*self.args,**self.kwargs)
        self.kwargs['fname'] = file_gen_new(fd.fname,fextend=fd.ftype)
        self.kwargs['energy'] = self.overall_energy
        fd = SaveFile(self.overall_system,*self.args,**self.kwargs)
        fd.run()
        outfile = fd.fname
        print('Note: file is saved to < {:} >'.format(outfile))

        filedict = {}

        # images
        if self.boim:
            filedict['probability data file'] = self.save_probdata()

            if len(self.overall_prob_begin['ball']):
                fgp = file_gen_new('bonds-all',fextend='png',foriginal=False)
                fbo = plot_save_image(
                    self.overall_prob_begin['ball'],
                    self.overall_prob_final['ball'],
                    dt=self.btol,
                    fname=fgp,
                    key='bonds',
                )
                if fbo: filedict['image all bonds filtration file'] = fgp

            if len(self.overall_prob_begin['aall']):
                fgp = file_gen_new('angles-all',fextend='png',foriginal=False)
                fbo = plot_save_image(
                    self.overall_prob_begin['aall'],
                    self.overall_prob_final['aall'],
                    dt=self.atol,
                    fname=fgp,
                    key='angles',
                )
                if fbo: filedict['image all angles filtration file'] = fgp

            filedict['image bonds par filtration file'] = []
            for i,t in enumerate(self.overall_prob_begin['bpar']):
                mark = 'bonds-par-{:}+{:}'.format(*self.bcon[i])
                fgp = file_gen_new(mark,fextend='png',foriginal=False)
                fbo = plot_save_image(
                    t,
                    self.overall_prob_final['bpar'][i],
                    dt=self.btol,
                    fname=fgp,
                    key='bonds',
                )
                if fbo: filedict['image bonds par filtration file'].append(fgp)

            filedict['image angles par filtration file'] = []
            for i,t in enumerate(self.overall_prob_begin['apar']):
                mark = 'angles-par-{:}+{:}+{:}'.format(*self.acon[i])
                fgp = file_gen_new(mark,fextend='png',foriginal=False)
                fbo = plot_save_image(
                    t,
                    self.overall_prob_final['apar'][i],
                    dt=self.btol,
                    fname=fgp,
                    key='angles',
                )
                if fbo: filedict['image angles par filtration file'].append(fgp)

        ftot = file_gen_new('bulk-process-info')
        print('Note: please check summary file for more info: < {:} >'.format(ftot))
        with open(ftot,'wt') as f:
            f.write('Note: starting  time: {:}\n'.format(time.ctime(self.mytime)))
            now = time.time()
            f.write('Note: finishing time: {:}\n'.format(time.ctime(now)))
            minutes = round((now-self.mytime)/60.0,2)
            f.write('Note: execution time: {:} minutes\n'.format(minutes))
            f.write('Note: current work path:\n')
            f.write('  => {:}\n'.format(os.path.abspath('.')))
            f.write('Note: random seed: {:}\n'.format(self.seed))
            f.write('Note: bulk process for input files:\n')
            for i,fd in enumerate(self.datafilelist):
                f.write('  => {:} -- molnms {:} => remove {:}\n'.format(fd,self.molnms[i],self.rmnmlist[i]))
            f.write('Note: total number of inputs: {:}\n'.format(sum(self.molnms)))
            if len(self.indexfilelist):
                f.write('\nNote: index files:\n')
                for fd in self.indexfilelist:
                    f.write('  => {:}\n'.format(fd))
            if self.mode is None or self.mode == 'dynamic':
                f.write('\nNote: filtration mode is: dynamic\n')
                if self.boall:
                    f.write('  => calculation is performed for all entries\n')
                else:
                    f.write('  => calculation is performed separately\n')
            else:
                f.write('\nNote: filtration mode is: static\n')
                if self.vndx: f.write('  => index value is: {:}\n'.format(self.vndx))
                if self.borandom:
                    f.write('  => randomly filtration')
                else:
                    f.write('  => lowest-bit filtration')
            f.write('\nNote: result file:\n')
            f.write('  => {:} -- molnms {:}\n'.format(outfile,len(self.overall_system)))
            f.write('\nNote: filtration ratio: {:}\n'.format(ratio))
            f.write('\nNote: index of connections start at 1\n')
            f.write('\nNote: bonds connections:\n')
            tot = ''
            out = '  => '
            for tmp in self.bcon:
                tmp = [i+1 for i in tmp]
                out += '{:}, '.format(tmp)
                if len(out) >= 80:
                    tot += out.rstrip() + '\n'
                    out = '  => '
            if out != '  => ': tot += out.rstrip() + '\n'
            tot += '\n'
            f.write(tot)

            f.write('Note: angles connections:\n')
            tot = ''
            out = '  => '
            for tmp in self.acon:
                tmp = [i+1 for i in tmp]
                out += '{:}, '.format(tmp)
                if len(out) >= 80:
                    tot += out.rstrip() + '\n'
                    out = '  => '
            if out != '  => ': tot += out.rstrip() + '\n'
            tot += '\n'
            f.write(tot)

            f.write('Note: btol   : {:} Angstrom\n'.format(self.btol))
            f.write('Note: atol   : {:} degree\n'.format(self.atol))

            for k,v in filedict.items():
                if isinstance(v,str):
                    f.write('Note: {:<40}:   {:}\n'.format(k,v))
                elif len(v):
                    f.write('Note: {:}:\n'.format(k))
                    for i,j in enumerate(v):
                        f.write('  ==> {:>3}: {:}\n'.format(i+1,j))

    def save_probdata(self):
        def gen_outputs(prob,bcon,acon,key):
            def fout(pdata):
                txt = ''
                out = '{:}\n'.format(pdata[1])
                for p in pdata[0]:
                    txt += '{:}  '.format(p)
                    if len(txt) >= 80:
                        out += txt.strip() + '\n'
                        txt = ''
                out += txt.strip() + '\n\n'
                return out

            k = key.upper()
            contents = f'@{k}  BALL\n' + fout(prob['ball']) if len(prob['ball']) else ''
            if len(prob['bpar']):
                for n,data in enumerate(prob['bpar']):
                    contents += '@{:}  BPAR   {:}  {:}\n'.format(k,*bcon[n])
                    contents += fout(data)
            if len(prob['aall']): contents += f'@{k}  AALL\n' + fout(prob['aall'])
            if len(prob['apar']):
                for n,data in enumerate(prob['apar']):
                    contents += '@{:}  APAR   {:}  {:}  {:}\n'.format(k,*acon[n])
                    contents += fout(data)
            return contents

        fdata = file_gen_new('bulk-probability-data')
        print('Note: probability data is saved to < {:} >'.format(fdata))
        with open(fdata,'wt') as f:
            f.write('# Note: filtration process probability data\n\n')
            for tmp in self.datafilelist:
                f.write('@FILE {:}\n'.format(tmp))
            f.write('@BTOL   {:}\n@ATOL   {:}\n\n\n'.format(self.btol, self.atol))
            mbcon = [[i+1 for i in j] for j in self.bcon]
            macon = [[i+1 for i in j] for j in self.acon]
            f.write(gen_outputs(self.overall_prob_begin,mbcon,macon,'begin'))
            f.write(gen_outputs(self.overall_prob_final,mbcon,macon,'final'))
        return fdata

    def get_datalist(self,filelist):
        """return 4D list"""
        datalist = []
        energylist = []
        for f in filelist:
            rf = ReadFile(f)
            if rf.nice:
                rf.run()
                print('Note: for file < {:} >, number of inputs < {:} >'.format(f,len(rf.system)))
            else:
                print(rf.info)
            datalist.append(rf.system)
            energylist.append(rf.energy)
        return datalist,energylist

    def get_connections(self,system):
        """
        Rule:
            now, we think user input is always on the first priority,
            and bcon & acon at the latter place

            connection is got in sequence:
            1) if bcon or acon is set
            2) if fragments is set
            3) default: bcon=AnglePerception.nconb, acon=AnglePerception.ncona
        """
        fn = AnglePerception(system)
        if not fn.nice:
            self.nice = False
            self.info = fn.info
            return
        fn.run()

        if 'userinputs' in self.kwargs and self.kwargs['userinputs'] is True:
            self.kwargs['userinputs'] = True
        else:
            self.kwargs['userinputs'] = False

        bog = False
        if 'fragments' in self.kwargs and self.kwargs['fragments'] is not None:
            bog = True
            if not len(self.kwargs['fragments']): self.kwargs['fragments'] = fn.fragments
            if self.kwargs['userinputs']:
                fg = [[j-1 for j in i] for i in self.kwargs['fragments']]
                if min([min(i) for i in fg]) < 0:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: atom index should start at 1'
                    return
                self.kwargs['fragments'] = fg
            else:
                fg = self.kwargs['fragments']
            if max([max(i) for i in fg]) >= len(system):
                self.nice = False
                self.info = 'Fatal: wrong defined: atom index exceeding: < {:} >'.format(len(system))
                return
            if len(set([j for i in fg for j in i])) != sum([len(i) for i in fg]):
                self.nice = False
                self.info = 'Fatal: wrong defined: fragments: has repeats'
                return
            fbcon, facon = self.func_calc_connections(self.kwargs['fragments'])
        else:
            self.kwargs['fragments'] = fn.fragments

        # take care of special case, when input is None
        if 'bcon' in self.kwargs and self.kwargs['bcon'] is not None:
            tmpbcon = self.kwargs['bcon']
            if isinstance(self.kwargs['bcon'],list):
                self.kwargs['bcon'] = self.check_user_input_connections(
                    self.kwargs['bcon'],
                    fn.fragments,
                    self.kwargs['userinputs']
                )
            elif isinstance(self.kwargs['bcon'],str):
                contmp = tmpbcon.lower()
                if 'no' in contmp or 'non' in contmp or 'none' in contmp:
                    self.kwargs['bcon'] = []
                elif hasattr(fn,contmp) and 'b' in contmp and 'con' in contmp:
                    self.kwargs['bcon'] = getattr(fn,contmp)
                else:
                    self.kwargs['bcon'] = False
            else:
                self.kwargs['bcon'] = False
            if self.kwargs['bcon'] is False:
                self.nice = False
                self.info += '\nFatal: wrong defined: bcon: {:}'.format(tmpbcon)
                return
        elif bog:
            self.kwargs['bcon'] = fbcon
        else:
            if fn.fragments is None or len(fn.fragments) == 1:
                self.kwargs['bcon'] = fn.nconb
            else:
                self.kwargs['bcon'] = fn.fnconb

        if 'acon' in self.kwargs and self.kwargs['acon'] is not None:
            if isinstance(self.kwargs['acon'],list):
                self.kwargs['acon'] = self.check_user_input_connections(
                    self.kwargs['acon'],
                    fn.fragments,
                    self.kwargs['userinputs']
                )
            elif isinstance(self.kwargs['acon'],str):
                tmpacon = self.kwargs['acon']
                contmp = tmpacon.lower()
                if 'no' in contmp or 'non' in contmp or 'none' in contmp:
                    self.kwargs['acon'] = []
                elif hasattr(fn,contmp) and 'a' in contmp and 'con' in contmp:
                    self.kwargs['acon'] = getattr(fn,contmp)
                else:
                    self.kwargs['acon'] = False
            else:
                tmpacon = self.kwargs['acon']
                self.kwargs['acon'] = False
            if self.kwargs['acon'] is False:
                self.nice = False
                self.info = 'Fatal: wrong defined: acon: {:}'.format(tmpacon)
                return
        elif bog:
            self.kwargs['acon'] = facon
        else:
            if fn.fragments is None or len(fn.fragments) == 1:
                self.kwargs['acon'] = fn.ncona
            else:
                self.kwargs['acon'] = fn.fncona

        # always set final userinputs to False
        self.kwargs['userinputs'] = False

    def check_user_input_connections(self,ul,fl,bo=None):
        if not len(ul): return []
        offset = 1 if bo is True else 0
        for ndx in ul:
            if len(set(ndx)) != len(ndx):
                self.nice = False
                self.info = 'Fatal: wrong defined: repeats at: < {:} >'.format(ndx)
                return False
        ux = max([max(i) for i in ul])
        fx = max([max(i) for i in fl])
        if ux-offset > fx:
            self.nice = False
            self.info = 'Fatal: wrong defined: atom index exceeding: < {:} >'.format(ux)
            return False
        # convert to python readable number, starts at 0
        ptmp = []
        for i in ul: ptmp.append([j-offset for j in i])
        for i in ptmp:
            if min(i) < 0:
                self.nice = False
                self.info = 'Fatal: atom index should start at 1'
                return False
        return ptmp

    def func_calc_connections(self,fragments):
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


class PlotSamples(BulkProcess):
    def __init__(self,probdatafilelist=None,nmsamples=None,nmlist=None,
                startndx=None,endndx=None,incndx=None,nmranges=None,
                *args,**kwargs):
        seed = kwargs['seed'] if 'seed' in kwargs else None
        self.seed = seed if seed else random.randrange(100000000)
        random.seed(self.seed)
        kwargs['seed'] = self.seed
        super().__init__(*args,**kwargs)
        if not len(self.datafilelist):
            self.nice = True
            self.info = ''

        self.probdatafilelist = []
        if probdatafilelist is not None:
            for f in probdatafilelist:
                if os.path.isfile(f):
                    self.probdatafilelist.append(f)
                else:
                    print('Warning: not probability data file < {:} >, ignoring'.format(f))

        if not len(self.datafilelist) and not len(self.probdatafilelist):
            self.nice = True
            self.info = 'Fatal: no valid inputs: datafilelist/probdatafilelist'
            return

        if 'userinputs' in kwargs and kwargs['userinputs']:
            if startndx is not None: startndx -= 1
            if endndx is not None: endndx -= 1
        
        # assume data has been properly processed, index starts at 0
        self.nmsamples = nmsamples
        self.startndx = startndx
        self.endndx = endndx
        self.incndx = incndx
        self.nmranges = nmranges
        self.nmlist = nmlist

        self.choices = []
        # after debug run, everything is ready
        allsystems = []
        if self.nice and len(self.datafilelist):
            allsystems = super().run(debug=True)
            if allsystems is None: allsystems = []
            tot = len(allsystems)
        if not len(self.datafilelist) or not self.nice:
            pass
        elif nmlist is None:
            pass
        elif isinstance(nmlist,list):
            for i in nmlist:
                if not isinstance(i,int):
                    self.nice = False
                    self.info = 'Fatal: wrong defined: not a number: {:}'.format(i)
                    return
                if i > tot:
                    self.nice = False
                    self.info = 'Fatal: wrong defined: too large: {:}'.format(i)
                    return
        else:
            self.nice = False
            if nmlist is None: self.info = 'Fatal: wrong defined: nmlist'
        # nmlist is in the highest priority
        bo = True
        if len(self.datafilelist) and self.nice and nmlist is not None:
            bo = False
            for i in nmlist:
                self.choices.append(random.sample(allsystems,i))
        if len(self.datafilelist) and bo and self.nice and 'datafilelist' in kwargs:
            if tot <= 20:
                self.info = 'Warning: too few inputs: datafilelist'
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if endndx is not None and endndx > tot:
                self.info = 'Fatal: too large: endndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if startndx is not None and startndx > tot:
                self.info = 'Fatal: too large: startndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if incndx is not None and incndx > tot:
                self.info = 'Fatal: too large: incndx --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if nmranges is not None and nmranges > tot:
                self.info = 'Fatal: too large: nmranges --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            if nmsamples is not None and nmsamples*8 > tot:
                self.info = 'Fatal: too large: nmsamples --> total:{:}'.format(tot)
                self.nice = False
        if len(self.datafilelist) and bo and self.nice:
            # alias
            samples = allsystems
            if startndx is not None: samples = samples[startndx:]
            if endndx is not None: samples = samples[:endndx]
            if nmranges is None: nmranges = 0
            if nmsamples is None: nmsamples = 5
            tot = len(samples)
            if incndx is None:
                t = tot // nmsamples
                for i in range(nmsamples):
                    self.choices.append(random.sample(samples,t))
            else:
                dt = random.randrange(nmranges)
                while dt < tot:
                    while True:
                        tmp = random.randrange(nmranges) * random.choice([-1,0,1])
                        if incndx + tmp > 0: break
                    if dt + incndx + tmp > 0: break
                    self.choices.append(samples[dt:dt+incndx+tmp])
                    dt += incndx + tmp
            # to avoid waste of time,
            # only number of samples bigger than 10 will be kept
            self.choices = [i for i in self.choices if len(i) > 10]

    def run(self):
        bondsdict = {'all':[], }
        anglesdict = {'all':[], }
        if len(self.probdatafilelist):
            print('Note: generate plots on probability data files ...')
            for f in self.probdatafilelist:
                print('Note: reading data from file: < {:} >'.format(f))
                bonds, angles = self.read_probdatafile(f)
                for k in bonds:
                    if k not in bondsdict: bondsdict[k] = []
                    bondsdict[k].append(bonds[k])
                for k in angles:
                    if k not in anglesdict: anglesdict[k] = []
                    anglesdict[k].append(angles[k])

        if len(self.choices):
            overall_prob_begin = []
            overall_prob_final = []
            overall_dt = []
            bconlist = []
            aconlist = []
            for cnt,ms in enumerate(self.choices):
                print('Note: generate plots on < {:} > sample files ...'.format(cnt+1))
                mf = Filtration(system=ms,*self.args,**self.kwargs)
                mf.run()
                bconlist.append(mf.bcon)
                aconlist.append(mf.acon)
                overall_prob_begin.append(mf.prob_begin)
                overall_prob_final.append(mf.prob_final)
                overall_dt.append([mf.btol, mf.atol])
            
            for i,di in enumerate(overall_prob_begin):
                df = overall_prob_final[i]

                if len(di['ball']):
                    p = {}
                    p['begin'] = di['ball']
                    p['final'] = df['ball']
                    p['dt'] = overall_dt[i][0]
                    bondsdict['all'].append(p)
                
                if len(di['aall']):
                    p = {}
                    p['begin'] = di['aall']
                    p['final'] = df['aall']
                    p['dt'] = overall_dt[i][1]
                    anglesdict['all'].append(p)
                
                if bconlist[i] is not None and not len(bconlist[i]) and not len(di['bpar']):
                    for j,k in enumerate(bconlist[i]):
                        if k[0] > k[1]: k[0],k[1] = k[1],k[0]
                        key = '{:}-{:}'.format(k[0],k[1])
                        if key not in bondsdict: bondsdict[key] = []
                        if len(di['bpar'][j]):
                            p = {}
                            p['begin'] = di['bpar'][j]
                            p['final'] = df['bpar'][j]
                            p['dt'] = overall_dt[i][0]
                            bondsdict[key].append(p)

                if aconlist[i] is not None and not len(aconlist[i]) and not len(di['apar']):
                    for j,k in enumerate(aconlist[i]):
                        if k[0] > k[2]: k[0],k[2] = k[2],k[0]
                        key = '{:}-{:}-{:}'.format(k[0],k[1],k[2])
                        if key not in anglesdict: anglesdict[key] = []
                        if len(di['apar'][j]):
                            p = {}
                            p['begin'] = di['apar'][j]
                            p['final'] = df['apar'][j]
                            p['dt'] = overall_dt[i][1]
                            anglesdict[key].append(p)
        # more info
        print('\nNote: random seed: {:}'.format(self.seed))
        if len(self.choices):
            if 'mode' in self.kwargs and self.kwargs['mode'] is not None:
                mode = self.kwargs['mode'].lower()
                mode = 'dynamic' if mode in ['d','dynamic'] else 'static'
            else:
                mode = 'dynamic'
            print('Note: filtration mode is: {:}'.format(mode))
            if mode == 'dynamic':
                boall = self.kwargs['boall'] if 'boall' in self.kwargs else None
                boall = True if boall or boall is None else False
                if boall:
                    print('  => calculation is performed for all entries')
                else:
                    print('  => calculation is performed separately')
            else:
                vndx = self.kwargs['vndx'] if 'vndx' in self.kwargs else None
                if vndx: print('  => index value is: {:}'.format(vndx))
                borandom = self.kwargs['borandom'] if 'borandom' in self.kwargs else None
                if borandom:
                    print('  => randomly filtration')
                else:
                    print('  => lowest-bit filtration')
        for k in bondsdict:
            self.save_image_samples(bondsdict[k],label='bonds')
        for k in anglesdict:
            self.save_image_samples(anglesdict[k],label='angles')

    class CYCLE:
        """infinite cycle works similar like itertools.cycle()"""
        def __init__(self,data):
            self.data = data
            self.nm = 0
        def __iter__(self):
            return self
        def __next__(self):
            if self.nm >= len(self.data): self.nm = 0
            v = self.data[self.nm]
            self.nm += 1
            return v
        def next(self):
            return self.__next__()

    def save_image_samples(self,datalist,label=None,fname=None):
        """
        Input:
            datalist: List[dict]: key in dict begin|final|dt
        """
        # the number of linestyle should be in ODD number
        linestyle = self.CYCLE(['solid','dotted','dashdot',])
        colors = self.CYCLE(['b','r','m','g','y','brown','palegreen','deepskyblue'])
        if label is None or label.lower() == 'bonds':
            label = 'bonds'
        else:
            label = 'angles'
        if fname is None: fname = 'bulk-image-samples-{:}'.format(label)
        fname = file_gen_new(fname,fextend='png')
        info = 'Conformers filtration on {:}'.format(label.capitalize())

        prodatalist = []
        molnms = []
        for data in datalist:
            if not isinstance(data,dict): continue
            if not ('begin' in data and 'final' in data and 'dt' in data): continue
            if len(data['begin'][0]) < 5 or len(data['final'][0]) < 5: continue
            prodatalist.append(data)
            molnms.append(sum(data['begin'][0]))
        if not len(prodatalist): return

        # sort them to make plot more tight
        reflist = sorted(range(len(molnms)), key=lambda k: molnms[k])
        prodatalist = [prodatalist[i] for i in reflist]

        # space [1,1] for plot, space [1,2] for legend
        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4,1]})
        fig.set_figheight(6)
        fig.set_figwidth(10)

        lines = []
        for d in prodatalist:
            color = colors.next()
            ils = linestyle.next()
            fls = linestyle.next()
            itxt = 'Begin molnms = {:}'.format(sum(d['begin'][0]))
            ftxt = 'Final molnms = {:}'.format(sum(d['final'][0]))
            ix = [d['begin'][1]+d['dt']*i for i in range(len(d['begin'][0]))]
            fx = [d['final'][1]+d['dt']*i for i in range(len(d['final'][0]))]
            iln, = ax1.plot(ix,d['begin'][0],color=color,linestyle=ils,label=itxt)
            fln, = ax1.plot(fx,d['final'][0],color=color,linestyle=fls,label=ftxt)
            lines.append(iln)
            lines.append(fln)

        ax1.set_title(info)
        ax2.axis('off')
        ax2.legend(handles=lines,loc='center')
        plt.tight_layout()
        print('Note: image file is saved to < {:} >'.format(fname))
        plt.savefig(fname)
        plt.close()
        return fname

    def read_probdatafile(self,file):
        def getdata(text,key=None):
            ltmp = text.split()
            if not len(ltmp): return []
            if key is not None and key >= len(ltmp): return []
            value = None
            data = []
            for i,v in enumerate(ltmp):
                if i == key:
                    try:
                        value = float(v)
                    except ValueError:
                        return []
                    continue
                else:
                    try:
                        v = int(v)
                    except ValueError:
                        return []
                data.append(v)
            return [data,value]

        atol = None
        btol = None
        ball = {'begin':'', 'final':''}
        bpar = {'begin':[], 'final':[]}
        aall = {'begin':'', 'final':''}
        apar = {'begin':[], 'final':[]}
        with open(file,'rt') as f: profile = f.readlines()
        cnt = -1
        while True:
            cnt += 1
            if cnt >= len(profile): break
            ltmp = profile[cnt].split()
            if len(ltmp) < 2 or ltmp[0][0] == '#': continue
            key = ltmp[0].lower()
            if key == '@btol' or key == '@atol':
                try:
                    if key == '@btol':
                        btol = float(ltmp[1])
                    else:
                        atol = float(ltmp[1])
                except ValueError:
                    continue
            elif key not in ['@begin','@final']:
                continue
            mark = ltmp[1].lower()
            if mark in ['ball','aall','bpar','apar']:
                txt = ''
                if mark == 'bpar':
                    if len(ltmp) < 4: continue
                    txt += ltmp[2] + '  ' + ltmp[3]
                if mark == 'apar':
                    if len(ltmp) < 5: continue
                    txt += ltmp[2] + '  ' + ltmp[3] + '  ' + ltmp[4]
                while True:
                    cnt += 1
                    if cnt >= len(profile): break
                    line = profile[cnt].strip()
                    if not len(line): continue
                    if line[0] == '@':
                        # caution: here is important
                        cnt -= 1
                        break
                    txt += '  ' + line
                if key == '@begin':
                    if mark == 'ball':
                        ball['begin'] = txt
                    elif mark == 'aall':
                        aall['begin'] = txt
                    elif mark == 'bpar':
                        bpar['begin'].append(txt)
                    else:
                        apar['begin'].append(txt)
                else:
                    if mark == 'ball':
                        ball['final'] = txt
                    elif mark == 'aall':
                        aall['final'] = txt
                    elif mark == 'bpar':
                        bpar['final'].append(txt)
                    else:
                        apar['final'].append(txt)
        if btol is None: btol = 0.1
        if atol is None: atol = 0.1
        # process for plot
        # bonds format:
        #   b1-b2 (b1<b2):   begin: List[ List[int], float]
        #                    final: List[ List[int], float]
        #                    dt   : float
        # angles format:
        #   a1-a-a2 (a1<a2): begin: List[ List[int], float]
        #                    final: List[ List[int], float]
        #                    dt   : float
        # special key:  'all'
        bonds = {'all':{}, }
        angles = {'all':{}, }

        bonds['all']['begin'] = getdata(ball['begin'],0)
        bonds['all']['final'] = getdata(ball['final'],0)
        bonds['all']['dt'] = btol
        angles['all']['begin'] = getdata(aall['begin'],0)
        angles['all']['final'] = getdata(aall['final'],0)
        angles['all']['dt'] = atol

        for i in bpar:
            for j in bpar[i]:
                data = getdata(j,2)
                if not len(data): continue
                bi,bj = data[0][0], data[0][1]
                if bi > bj: bj, bi = bi, bj
                key = '{:}-{:}'.format(bi,bj)
                data[0] = data[0][2:]
                if key not in bonds: bonds[key] = {}
                bonds[key][i] = data
                bonds[key]['dt'] = btol

        for i in apar:
            for j in apar[i]:
                data = getdata(j,3)
                if not len(data): continue
                ai,aj = data[0][0], data[0][2]
                if ai > aj: aj, ai = ai, aj
                key = '{:}-{:}-{:}'.format(ai,data[0][1],aj)
                data[0] = data[0][3:]
                if key not in angles: angles[key] = {}
                angles[key][i] = data
                angles[key]['dt'] = atol

        return bonds, angles


def parsecmd():
    """Parse command line input"""
    def parse_remove_chars(line):
        line = line.replace('[',' ').replace(']',' ').replace(';',',')
        return line

    def parse_in_line(line,bobcon=False,boacon=False):
        line = parse_remove_chars(line)
        line = line.replace('-',' ').replace('_',' ')
        lt = line.split(',')
        ls = []
        for t in lt:
            lp = t.split()
            if not len(lp): continue
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

        if not len(ls): return []

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
        '-f','--datafilelist',
        help='Data files, separate by space or comma',
        nargs='+',
        metavar='file',
    )
    parser.add_argument(
        '-d','--indexfilelist',
        help='Index files, as the filtration reference',
        nargs='+',
        metavar='file',
    )
    parser.add_argument(
        '--static',
        help='turn on static mode calculation, default is in dynamic/all mode',
        action='store_true',
    )
    parser.add_argument(
        '--separate',
        help='valid in dynamic mode, change to dynamic/separate mode',
        action='store_true',
    )
    parser.add_argument(
        '--vndx',
        help='valid in static mode, set filtration index value',
        type=float,
    )
    parser.add_argument(
        '--borandom',
        help='valid in static mode, rather than by lowest-bit, filtering out molecules randomly',
        action='store_true',
    )
    parser.add_argument(
        '--seed',
        help='random seed, (optional, highest priority)',
        type=int,
    )
    parser.add_argument(
        '-btol','--btol',
        help='bonds tolerance, any changes in bcon smaller than it will be filtered out',
        type=float,
    )
    parser.add_argument(
        '-atol','--atol',
        help='angles tolerance, any changes in acon smaller than it will be filtered out',
        type=float,
    )
    parser.add_argument(
        '-bcon',
        help='Bond connections, in pairs, separate by comma',
        nargs='+',
        metavar='B1 B2, B1-B2',
    )
    parser.add_argument(
        '-acon',
        help='Angle connections, in pairs, separate by comma',
        nargs='+',
        metavar='A1 A2 A3, A1-A2-A3',
    )
    parser.add_argument(
        '-g','--fragments',
        help='fragments, separate by comma',
        nargs='+',
        metavar='G',
    )
    parser.add_argument(
        '--no-oball',
        help='turn off bonds probability overall calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-obpar','--obpar',
        help='turn on bonds probability parameters calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '--no-oaall',
        help='turn off angles probability overall calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-oapar','--oapar',
        help='turn on angles probability parameters calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-nc','--no-force-double-check',
        help='turn off double check prompt info before execution',
        action='store_true',
    )
    parser.add_argument(
        '--features',
        help='show development features',
        action='store_true',
    )
    parser.add_argument(
        '-p','--file-format-explanations',
        help='show input system file format explanations',
        action='store_true',
    )
    parser.add_argument(
        '-nu','--no-userinputs',
        help='specify the indexes of input connections start at 0',
        action='store_true',
    )
    parser.add_argument(
        '-o','--fname',
        help='output system file name',
    )
    parser.add_argument(
        '-ft','--ftype',
        help='output system file type, [txt, xsf, xyz]',
    )
    subparser = parser.add_subparsers(title='continuous subcommand')
    sub = subparser.add_parser(
        'plot',
        help='plot cross comparsions based on different chosen samples',
        allow_abbrev=False
    )
    sub.set_defaults(command='plot')
    sub.add_argument(
        '-bf','--probdatafilelist',
        help='bulk process probability files',
        metavar='B',
        nargs='+',
    )
    sub.add_argument(
        '-nl','--nmlist',
        help='highest priority, select number of samples for plot, (optional)',
        metavar='n',
        nargs='+',
        type=int,
    )
    sub.add_argument(
        '-ns','--nmsamples',
        help='choose number of samples for plot, default is 5',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-sn','--startndx',
        help='start index for the inputs, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-en','--endndx',
        help='end index for the inputs, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-inc','--incndx',
        help='increments for choose, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-nr','--nmranges',
        help='random ranges for increments, (optional)',
        metavar='n',
        type=int,
    )

    # annoying part
    bo = False
    if len(sys.argv) == 1: bo = True
    if not bo and 'plot' not in sys.argv:
        if '-h' in sys.argv or '--help' in sys.argv: bo = True
    if bo:
        parser.print_help()
        exit()

    # due to the potential bug, args may not have enough namespace
    args,left = sub.parse_known_args()
    if 'plot' in left:
        left.remove('plot')
        opts = parser.parse_args(left)
        # conclude them
        args = argparse.Namespace(**vars(opts),**vars(args))
    else:
        args = parser.parse_args(sys.argv[1:])

    if 'features' in args and args.features:
        for i in FEATURES:
            print(i)
        exit()

    if 'file_format_explanations' in args and args.file_format_explanations:
        print(FILEFORMAT)
        exit()

    # default settings
    fdict = {
        'datafilelist'              :   None,
        'indexfilelist'             :   None,
        'bcon'                      :   None,
        'btol'                      :   None,
        'oball'                     :   True,
        'obpar'                     :   False,
        'acon'                      :   None,
        'atol'                      :   None,
        'oaall'                     :   True,
        'oapar'                     :   False,
        'mode'                      :   None,
        'vndx'                      :   None,
        'borandom'                  :   None,
        'boall'                     :   True,
        'fragments'                 :   None,
        'bool_force_double_check'   :   True,
        'userinputs'                :   True,
        'fname'                     :   None,
        'ftype'                     :   None,
        'probdatafilelist'          :   None,
        'nmlist'                    :   None,
        'nmsamples'                 :   None,
        'startndx'                  :   None,
        'endndx'                    :   None,
        'incndx'                    :   None,
        'nmranges'                  :   None,
        'seed'                      :   None,
    }

    bod = False
    if 'datafilelist' in args and args.datafilelist is not None:
        if len(args.datafilelist): bod = True
    bog = False
    if 'probdatafilelist' in args and args.probdatafilelist is not None:
        if len(args.probdatafilelist): bog = True
    if 'command' in args:
        if not bod and not bog:
            print('Fatal: no input: missing:  -f/--datafilelist or -bf/--probdatafilelist')
            exit()
    else:
        if not bod:
            print('Warning: -f/--datafilelist is missing')
            exit()
    
    if bod:
        stmp = parse_remove_chars(' '.join(args.datafilelist)).replace(',',' ')
        fdict['datafilelist'] = stmp.split()
    
    if bog:
        stmp = parse_remove_chars(' '.join(args.probdatafilelist)).replace(',',' ')
        fdict['probdatafilelist'] = stmp.split()

    if 'indexfilelist' in args and args.indexfilelist:
        stmp = parse_remove_chars(' '.join(args.indexfilelist)).replace(',',' ')
        fdict['indexfilelist'] = stmp.split()

    if 'static' in args and args.static: fdict['mode'] = 'static'
    if 'separate' in args and args.separate: fdict['boall'] = False
    if 'vndx' in args: fdict['vndx'] = args.vndx     # be aware of 0.0
    if 'borandom' in args and args.borandom: fdict['borandom'] = True

    if 'bcon' in args and args.bcon:
        fdict['bcon'] = parse_in_line(' '.join(args.bcon),bobcon=True)
    if 'acon' in args and args.acon:
        fdict['acon'] = parse_in_line(' '.join(args.acon),boacon=True)
    if 'fragments' in args and args.fragments:
        fdict['fragments'] = parse_in_line(' '.join(args.fragments))

    if 'btol' in args and args.btol: fdict['btol'] = args.btol
    if 'atol' in args and args.atol: fdict['atol'] = args.atol
    if 'oball' in args and args.no_oball: fdict['oball'] = False
    if 'obpar' in args and args.obpar: fdict['obpar'] = True
    if 'oaall' in args and args.no_oaall: fdict['oaall'] = False
    if 'oapar' in args and args.oapar: fdict['oapar'] = True
    if 'no_force_double_check' in args and args.no_force_double_check:
        fdict['bool_force_double_check'] = False
    if 'no_userinputs' in args and args.no_userinputs: fdict['userinputs'] = False
    if 'fname' in args and args.fname: fdict['fname'] = args.fname
    if 'ftype' in args and args.ftype: fdict['ftype'] = args.ftype
    if 'nmlist' in args and args.nmlist: fdict['nmlist'] = args.nmlist
    if 'nmsamples' in args and args.nmsamples: fdict['nmsamples'] = args.nmsamples
    if 'startndx' in args and args.startndx: fdict['startndx'] = args.startndx
    if 'endndx' in args and args.endndx: fdict['endndx'] = args.endndx
    if 'incndx' in args and args.incndx: fdict['incndx'] = args.incndx
    if 'nmranges' in args and args.nmranges: fdict['nmranges'] = args.nmranges
    if 'seed' in args and args.seed: fdict['seed'] = args.seed

    print('Note: time: {:}'.format(time.ctime()))
    if 'command' in args:
        print('Note: processing plot ...')
        PS = PlotSamples(**fdict)
    else:
        print('Note: processing data files ...')
        PS = BulkProcess(**fdict)
    if not PS.nice:
        print(PS.info)
        return
    PS.run()
    if not PS.nice: print(PS.info)


if __name__ == '__main__':
    parsecmd()


