#!/usr/bin/env python3

import os
import sys
import math
import argparse


FEATURES = [
    'version 0.1    : start',
    'version 0.2    : add support for save to com file',
    'version 0.2.2  : temporary, exclusive for Henry reaction',
    'version 0.2.3  : temporary, exclusive for TS SN2 reaction',
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


class ReadFile:
    """
    Initialization:
        file :
            format:
                txt file
                xsf file
                xyz file

    Attributes:
        system : System[Mol[Atom, ...]] : List[List[[atomtype, x,y,z], ...]]

        energy : 1D List[float]  :   None means not exist
    """
    def __init__(self,file,ext=None):
        self.file = file

        # decide file format
        if ext is None:
            ndx = file.rfind('.')
            if ndx == -1 or ndx + 1 >= len(file):
                self.ext = None
            else:
                self.ext = file[ndx+1:].lower()
            if self.ext is None: self.ext = 'txt'
        else:
            self.ext = ext
        
        if self.ext not in ['txt','xsf','xyz']:
            print('Input file format: {:}'.format(ext))
            raise ValueError('Error: input file format is not supported')



    def run(self):
        print('Note: reading file < {:} > ...'.format(self.file))
        self.energy = []
        if self.ext == 'txt':
            self.read_txt()
        elif self.ext == 'xsf':
            self.read_xsf()
        elif self.ext == 'xyz':
            self.read_xyz()
        else:
            print('UPDATENEEDED')
            self.system = []



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
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                prolist.append(ls)
                enelist.append(ene)

        self.system,self.energy = self.check_atomtype(prolist,enelist)



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
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                prolist.append(ls)
                enelist.append(ene)

        self.system,self.energy = self.check_atomtype(prolist,enelist)



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
                print('\nWarning: line {:}: {:}'.format(errnum,errline))
                print('Note: whole sets are omitted...\n')
            else:
                prolist.append(ls)
                enelist.append(ene)

        self.system,self.energy = self.check_atomtype(prolist,enelist)



    def check_atomtype(self,syslist,enelist):
        """check atomtype for syslist & return results basis on first entry

        Inputs:
            syslist : System[Mol[Atoms]] : List[List[[atomtype, x, y, z], ...]]
            enelist : List[Float or None]

        Return:
            same format with inputs
        """
        if len(syslist) <= 1: return syslist,enelist

        ndxlist = [i[0] for i in syslist[0]]
        prosyslist = [syslist[0], ]
        proenelist = [enelist[0], ]
        lth = len(ndxlist)
        for n,mol in enumerate(syslist[1:]):
            if len(mol) != lth:
                print('Wrong: not cooresponded: {:}'.format(mol[0]))
                print('Note: whole sets are omitted...\n')
                continue
            bo = True
            for cnt,atom in enumerate(mol):
                if ndxlist[cnt] != atom[0]:
                    print('Wrong: not cooresponded: {:}'.format(atom))
                    print('Note: whole sets are omitted...\n')
                    bo = False
                    break
                if not bo: break
            if bo:
                prosyslist.append(mol)
                # be careful with n starts
                proenelist.append(enelist[n+1])
        
        # make energy label more compatible
        totene = []
        for e in proenelist:
            if e is not None and abs(e) <= 0.5: e = None
            totene.append(e)

        return prosyslist,totene




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




class SaveFile:
    """opposite operation to ReadFile

    Inputs:
        system : 2D  :  Mol[Atom, ...]  :  List[[atomtype, x,y,z], ...]

        ftype  : Output file type
            format
                txt file
                xsf file
                xyz file
                com file

        force_real_atom_type    :   Boolean     :   default False
        force_double_check      :   Boolean     :   default True
        outfile_with_energy     :   True, False, None(all)

        fname : file to be saved, warning, overwritten may happen
    """
    def __init__(self,system,*args,**kwargs):
        self.system = system

        self.basis_set = '# HF/6-31G Pop=CHelpG'
        if 'basis_set' in kwargs and kwargs['basis_set'] is not None:
            self.basis_set = kwargs['basis_set']
        self.charge_spin = '0 1'
        if 'charge_spin' in kwargs and kwargs['charge_spin'] is not None:
            self.charge_spin = kwargs['charge_spin']

        # Rule
        # self.ftype always takes the precedence
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
            print('Warning: current not support < {:} >'.format(self.ftype))
            raise ValueError('not support')

        if self.fname is None: self.fname = 'system'
        # keep dot conversion in original fname, if it has
        ndx = self.fname.rfind('.')
        if ndx == -1:
            self.fname = self.fname + '.' + self.ftype
        else:
            ext = self.fname[ndx:]
            if ext != '.' + self.ftype:
                self.fname = self.fname + '.' + self.ftype
        
        self.fname = self.fname[:self.fname.rfind('.')]

        ft = False
        if 'force_real_atom_type' in kwargs:
            ft = kwargs['force_real_atom_type']
        self.force_real_atom_type = True if ft is True else False

        fc = True
        if 'force_double_check' in kwargs:
            fc = kwargs['force_double_check']
        self.force_double_check = False if fc is False else True

        fo = None
        if 'outfile_with_energy' in kwargs:
            fo = kwargs['outfile_with_energy']
        if fo is True:
            self.outfile_with_energy = True
        elif fo is False:
            self.outfile_with_energy = False
        else:
            self.outfile_with_energy = None

        if 'energy' in kwargs and kwargs['energy'] is not None:
            self.energy = kwargs['energy']
        else:
            self.energy = []
        
        self.outfiles = []

        self.henry = False
        if 'henry' in kwargs and kwargs['henry'] is True:
            self.henry = True
            self.ftype = 'com'
        
        self.sn2 = False
        if 'sn2' in kwargs and kwargs['sn2'] is True:
            self.sn2 = True
            self.ftype = 'com'



    def run(self):
        """save to file"""
        if len(self.system) == 0:
            print('Warning: no inputs')
            return

        if self.force_double_check:
            if len(self.energy) == 0:
                totene = 0
            else:
                totene = sum([1 for i in self.energy if i is not None])
            totsys = len(self.system)
            if self.outfile_with_energy is True:
                fnum = totene
            elif self.outfile_with_energy is False:
                fnum = totsys - totene
            else:
                fnum = totsys

            print('\nCheck: total number of inputs: < {:} >'.format(totsys))
            print('Check: number of molecules have energy: < {:} >'.format(totene))
            atp = 'no' if self.force_real_atom_type is False else 'yes'
            print('Check: whether force to use real atom type: < {:} >'.format(atp))
            enp = 'no' if self.outfile_with_energy is False else 'yes'
            print('Check: whether output files contain energy < {:} >'.format(enp))
            print('Check: output file type: < {:} >'.format(self.ftype))
            print('Check: output file base name: < {:} >'.format(self.fname))
            print('Check: number of files will be generated: < {:} >'.format(fnum))
            print('Check: one sample will be:\n')

        atlist = self.get_atomtype()
        if self.ftype == 'txt':
            self.save_txt(atlist)
        elif self.ftype == 'xsf':
            self.save_xsf(atlist)
        elif self.ftype == 'xyz':
            self.save_xyz(atlist)
        elif self.ftype == 'com':
            self.save_com(atlist)
        else:
            # this will never be executed
            print('UPDATENEEDED')
        
        if len(self.outfiles) == 0:
            print('Warning: no data after processed')
            return
        if self.force_double_check:
            print(self.outfiles[0])
            print('Check: do you want to continue? y/yes, else not:',end='    ')
            tmp = input()
            if tmp.lower() not in ['y','yes']:
                print('\nNote: you have decided to quit...')
                return
            print('\nNote: files generating: < {:} >'.format(fnum))
        self.save_file()



    def save_file(self):
        # special case for Gaussian com files
        if self.ftype == 'com':
            for mol in self.outfiles:
                fout = file_gen_new(self.fname,fextend=self.ftype,foriginal=False)
                chk = '%chk={:}.chk'.format(fout[:fout.rfind('.')])
                mol = chk + mol[mol.find('\n'):]
                with open(fout,'wt') as f:
                    f.write(mol)
                    f.write('\n\n')
            return

        for mol in self.outfiles:
            fout = file_gen_new(self.fname,fextend=self.ftype,foriginal=False)
            with open(fout,'wt') as f:
                f.write(mol)
                f.write('\n\n')
        print('Note: done on files generation')



    def save_xsf(self,atlist):
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    fout = '#  {:}\n\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    fout = '#\n\n'
            if bo:
                fout += 'ATOMS\n'
                for cnt,at in enumerate(mol):
                    line = '{:<2} {:>12} {:>12} {:>12}   1.0  1.0  1.0\n'.format(atlist[cnt],*at[1:])
                    fout += line
                self.outfiles.append(fout)



    def save_txt(self,atlist):
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    fout = '#  {:}\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    fout = '#\n\n'
            if bo:
                for cnt,at in enumerate(mol):
                    fout += '{:<2} {:>12} {:>12} {:>12}\n'.format(atlist[cnt],*at[1:])
                self.outfiles.append(fout)



    def save_xyz(self,atlist):
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    line = 'Properties=species:S:1:pos:R:3 energy={:}\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    line = 'Properties=species:S:1:pos:R:3 energy=0.0\n'
            if bo:
                fout = '{:}\n'.format(len(mol))
                fout += line
                for cnt,at in enumerate(mol):
                    fout += '{:<2} {:>12} {:>12} {:>12}\n'.format(atlist[cnt],*at[1:])
                self.outfiles.append(fout)



    def save_com(self,atlist):
        if self.henry:
            self.save_com_henry(atlist)
            return
        if self.sn2:
            self.save_com_sn2(atlist)
            return
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    title = 'Warning: already calc energy={:}\n\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    title = 'BOSS Backupfile-* Process\n\n'
            if bo:
                fout = '%chk=file-{:}.chk\n'.format(ndx+1)
                fout += self.basis_set + '\n\n' + title + self.charge_spin + '\n'
                for cnt,at in enumerate(mol):
                    fout += '{:<2} {:>12} {:>12} {:>12}\n'.format(atlist[cnt],*at[1:])
                self.outfiles.append(fout)



    def save_com_henry(self,atlist):
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    title = 'Warning: already calc energy={:}\n\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    title = 'BOSS Backupfile-* Process\n\n'
            if bo:
                fout = '%chk=file-{:}.chk\n%NProcShared=1\n%Mem=4000MB\n'.format(ndx+1)
                fout += '# MP2/6-31+G(d) pop=CM5\n# SCF=InCore Transformation=InCore\n'
                fout += '\nHenry Reaction\n\n-1 1\n'
                for cnt,at in enumerate(mol):
                    fout += '{:<2} {:>12} {:>12} {:>12}\n'.format(atlist[cnt],*at[1:])
                self.outfiles.append(fout)



    def save_com_sn2(self,atlist):
        for ndx,mol in enumerate(self.system):
            bo = False
            if len(self.energy) != 0 and self.energy[ndx] is not None:
                if self.outfile_with_energy is not False:
                    bo = True
                    title = 'Warning: already calc energy={:}\n\n'.format(self.energy[ndx])
            else:
                if self.outfile_with_energy is not True:
                    bo = True
                    title = 'BOSS Backupfile-* Process\n\n'
            if bo:
                fout = '%chk=file-{:}.chk\n'.format(ndx+1)
                fout += '# M062X/6-31+G(d) Pop=CM5\n# SCF=InCore Transformation=InCore\n'
                fout += '\nTS SN2 Chlorine MethylChloride\n\n-1 1\n'
                for cnt,at in enumerate(mol):
                    fout += '{:<2} {:>12} {:>12} {:>12}\n'.format(atlist[cnt],*at[1:])
                self.outfiles.append(fout)



    def get_atomtype(self):
        pt_nm = {
            1:'H', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 15:'P', 16:'S', 17:'Cl',
            35:'Br', 53:'I',
        }
        pt_ch = {}
        for k,v in pt_nm.items(): pt_ch[str(k)] = v
        atlist = [j[0] for i in self.system for j in i]
        if self.force_real_atom_type:
            reflist = []
            for at in atlist:
                if at in pt_nm:
                    at = pt_nm[at]
                elif at in pt_ch:
                    at = pt_ch[at]
                reflist.append(at)
            atlist = reflist
        return atlist




def parsecmd():
    """Parse command line input"""
    parser = argparse.ArgumentParser(
        description='BOSS Backup-* File Process',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION,
    )
    parser.add_argument(
        '-f','--backupfile',
        help='BOSS Backup-* file',
    )
    parser.add_argument(
        '-nc','--no-force-double-check',
        help='Turn off double check before execution',
        action='store_true',
    )
    parser.add_argument(
        '-nt','--no-force-real-atom-type',
        help='No using real atom type, default is True, if possible',
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
        help='Output file base name',
    )
    parser.add_argument(
        '-ft','--ftype',
        help='Output file type, [txt, xsf, xyz, com]',
    )
    parser.add_argument(
        '-bs','--basis_set',
        help='Valid only on com file, to be compatible on SHELL parse, ' \
            'input pound sign (#) is @: e.g. "@ HF/6-31G" means "# HF/6-31G". ' \
            'Note: input basis_set & charge_spin will not be again validated, ' \
            'please take care. Default: "# HF/6-31G Pop=CHelpG"',
        nargs='+',
    )
    parser.add_argument(
        '-cs','--charge_spin',
        help='Valid only on com file, default: "0 1"',
        nargs='+',
    )
    parser.add_argument(
        '-oe','--outfile-with-energy',
        help='Filter, whether output files contain energy, [all, y/yes/true, n/no/false]'
    )
    parser.add_argument(
        '--henry',
        action='store_true',
        help='temporary fix, for Henry only'
    )
    parser.add_argument(
        '--sn2',
        action='store_true',
        help='temporary fix, for SN2 Chlorine MethylChloride only'
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
        'backupfile'            :   None,
        'force_double_check'    :   True,
        'force_real_atom_type'  :   True,
        'fname'                 :   None,
        'ftype'                 :   'xsf',
        'outfile_with_energy'   :   None,
        'basis_set'             :   None,
        'charge_spin'           :   None,
        'henry'                 :   False,
        'sn2'                   :   False,
    }
    if args.backupfile is None:
        print('Warning: -f/--backupfile is missing')
        exit()

    fdict['backupfile'] = args.backupfile
    if args.no_force_double_check: fdict['force_double_check'] = False
    if args.no_force_real_atom_type: fdict['force_real_atom_type'] = False
    if args.fname is not None: fdict['fname'] = args.fname
    if args.ftype is not None: fdict['ftype'] = args.ftype
    if args.henry: fdict['henry'] = True
    if args.sn2: fdict['sn2'] = True
    if args.outfile_with_energy is not None:
        if args.outfile_with_energy.lower() in ['y','yes','true']:
            fdict['outfile_with_energy'] = True
        elif args.outfile_with_energy.lower() in ['n','no','false']:
            fdict['outfile_with_energy'] = False
    
    if args.basis_set:
        bs = ' '.join(args.basis_set).replace('"', ' ').replace("'",' ')
        bs = bs.replace('@','#')
        fdict['basis_set'] = ' '.join(bs.split())
    if args.charge_spin:
        cs = ' '.join(args.charge_spin).replace('"', ' ').replace("'",' ')
        fdict['charge_spin'] = ' '.join(cs.split())
    
    
    RF = ReadFile(fdict['backupfile'])
    RF.run()

    fdict['energy'] = RF.energy
    SF = SaveFile(RF.system,**fdict)
    SF.run()




def read_txt(file):
    with open(file,mode='rt') as f:
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
    print(len(promol))



if __name__ == '__main__':
    parsecmd()


