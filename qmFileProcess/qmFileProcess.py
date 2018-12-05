import os
import argparse

def file_size_check(path,fsize=1000):
    """This function is used to check the file existence and size,
       the unit of size is in megabety"""
        
    try:
        sizetmp = os.stat(path).st_size
        if sizetmp/1024/1024 > fsize:
            print('Error: the file size is far larger than %f MB' % fsize)
            print('Error path: ',path)
            exit()
          
    except IOError:
        print('Error : cannot open the file!')
        print('Error path: ',path)
        exit()
        
    return 1



def file_gen_new(fname,fextend=None):
    """This function is used to make new file but without overriding the
       old one"""
    
    filename = fname
    pos = filename.rfind('.')
    if pos != -1:
        fname = filename[:pos]
        if fextend is None:
            fextend = filename[pos:]
        else:
            fextend = '.' + fextend
    else:
        if fextend is None:
            fextend = '.xsf'
        else:
            fextend = '.' + fextend
    
    i = 1
    filename = fname
    while True:
        fname = filename + '_' + str(i) + fextend
        try:
            f = open(fname)
            f.close()
        except:
            break
        i += 1
    return fname


class QMFile(object):
    def __init__(self,file,atomnm,fextend=None,fname=None,basis_set=None,gcheck=None,
                 charge_spin=None,gtitle=None):
        
        dump_value = file_size_check(file)        
        self.file = file

        if fname is None:
            self.fname = 'QM'
        else:
            self.fname = fname

        if fextend is None:
            pos = self.fname.rfind('.')
            if pos == -1:
                fextend = 'xsf'
            else:
                if pos + 1 == len(self.fname):
                    fextend = 'xsf'
                else:
                    fextend = self.fname[pos+1:]

        if fextend.lower() == 'xsf':
            self.fextend = 'xsf'
        elif fextend.lower() == 'pdb':
            self.fextend = 'pdb'
        elif fextend.lower() == 'com':
            self.fextend = 'com'
        else:
            print('Error: currently, only file formats `xsf, pdb, com\' are supported')
            exit()
            
        try:
            self.atomnm = int(atomnm)
            if self.atomnm <= 0:
                raise ValueError
        except:
            print('Error: the parameter atomnm has to be a positive integer number')
            exit()

        if basis_set is None:
            self.basis_set = '# HF/6-31G(d)'
        elif len(basis_set.split()) == 0:
            self.basis_set = '# HF/6-31G(d)'
        else:
            self.basis_set = basis_set

        if gcheck is None:
            self.gcheck = None
        elif len(gcheck.split()) == 0:
            self.gcheck = None
        else:
            self.gcheck = gcheck

        if charge_spin is None:
            self.charge_spin = '0 1'
        elif len(charge_spin.split()) == 0:
            self.charge_spin = '0 1'
        else:
            if len(charge_spin.split()) != 2:
                print('Error: wrong definition for the parameter charge_spin')
                print(charge_spin)
                exit()
            self.charge_spin = charge_spin

        if gtitle is None:
            self.gtitle = 'QM process'
        elif len(gtitle.split()) == 0:
            self.gtitle = 'QM process'
        else:
            self.gtitle = gtitle
            

        self.prolist = self.pro_qmfile()
        if len(self.prolist) == 0:
            print('Error: no inputs')
            exit()


    def file_print(self):
        periodic_table = { 1:'H', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F',
                           11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl',
                           19:'K', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn',
                           26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 35:'Br',
                           47:'Ag', 53:'I', 79:'Au', 80:'Hg', }
    
        atomlist = []
        for i in self.prolist[0]:
            atomlist.append(periodic_table[i[0]])

        fnamelist = []
        for mol in self.prolist:
            fname = file_gen_new(self.fname,self.fextend)
            fnamelist.append(fname)
            with open(fname,'wt') as f:
                if self.fextend == 'xsf':
                    f.write('#\n')
                    f.write('\n')
                    f.write('ATOMS\n')
                    ndx = 0
                    for i in mol:
                        line = '{:3} {:>10} {:>10} {:>10}   1.0  1.0  1.0\n'.format(atomlist[ndx],i[1],i[2],i[3])
                        f.write(line)
                        ndx += 1
                    f.write('\n')
                elif self.fextend == 'pdb':
                    f.write('REMARK   1 File created by QMFile process program\n')
                    ndx = 0
                    for i in mol:
                        line = 'ATOM  {:>5} {:>2}                {:>8}{:>8}{:>8}\n'.format(ndx+1,atomlist[ndx],i[1],i[2],i[3])
                        f.write(line)
                        ndx += 1
                    f.write('END\n\n')
                elif self.fextend == 'com':
                    if self.gcheck is None:
                        fcheck = fname[:fname.rfind('.')] + '.chk'
                        f.write('%chk={:}\n'.format(fcheck))
                    else:
                        f.write('%chk={:}\n'.format(self.gcheck))
                    f.write(self.basis_set)
                    f.write('\n\n')
                    f.write(self.gtitle)
                    f.write('\n\n')
                    f.write(self.charge_spin)
                    f.write('\n')
                    ndx = 0
                    for i in mol:
                        line = '{:4} {:>10} {:>10} {:>10}\n'.format(atomlist[ndx],i[1],i[2],i[3])
                        f.write(line)
                        ndx += 1
                    f.write('\n')
                        
                    
        pos = self.fname.rfind('.')
        filename = self.fname[:pos] if pos != -1 else self.fname
        
        pf = filename + '_NameList.txt'
        fname = file_gen_new(pf)
        with open(fname,'wt') as f:
            f.write('# This is a collection of all the generated files\n')
            f.write('#\n')
            f.write('# The input file is:\n')
            f.write('#   {:}\n'.format(self.file))
            f.write('# The total atoms number are: < {:} >\n\n\n'.format(self.atomnm))
            f.write('# The output files are;\n')
            for i in fnamelist:
                f.write(i)
                f.write('\n')
                

    def pro_qmfile(self):
        profile = []
        with open(self.file,mode='rt') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break
                else:
                    ltmp = line.split()
                    if len(ltmp) == 0:
                        continue
                    elif len(ltmp) == 4:
                        try:
                            atomid = int(ltmp[0])
                            rx = float(ltmp[1])
                            ry = float(ltmp[2])
                            rz = float(ltmp[3])
                        except:
                            print('Error: the input line is not correctly formated')
                            print(line)
                            exit()
                    else:
                        print('Error: the input line is wrong')
                        print(line)
                        exit()

                    profile.append([atomid,rx,ry,rz])
         
        if len(profile) % self.atomnm != 0:
            print('Error: the input file and atomnm are not corresponded')
            exit()

        # split to different molecules
        prolist = []
        i = 0
        while i < len(profile):
            j = 0
            ls =[]
            while j < self.atomnm:
                ls.append(profile[i+j])
                j += 1
            prolist.append(ls)
            i += self.atomnm

        # remove the repeats
        reflist = []
        i = 1
        while i < len(prolist):
            j = 0
            while j < i:
                bool_repeats = True
                for k in range(self.atomnm):
                    if prolist[i][k][0] != prolist[j][k][0]:
                        print('Error: the input file and atomnm are not corresponded')
                        exit()
                    for nm in [1,2,3]:
                        if prolist[i][k][nm] != prolist[j][k][nm]:
                            bool_repeats = False
                            break
                    if not bool_repeats:
                        break
                if bool_repeats:
                    reflist.append(i)
                    break
                j += 1
            i += 1


        ndxlist = list(range(len(prolist)))
        for i in reflist:
            j = ndxlist.index(i)
            prolist.remove(prolist[j])
            ndxlist.remove(i)

        return prolist

        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MC/FEP QM backup file process',allow_abbrev=False)
    parser.add_argument('-f','--file',help='Input file',required=True)
    parser.add_argument('-nm','--atomnm',help='Atom number',required=True)
    parser.add_argument('-o','--fname',help='Output file name, default is QM')
    parser.add_argument('-t','--fextend',help='Output file extension, this will be used as the final file extension,\
                        only xsf, pdb and com files are supported, default is xsf')
    parser.add_argument('-bs','--basis_set',help='Only Valid for `com\' file, default is `# HF/6-31G(d)\'',nargs='+')
    parser.add_argument('-cs','--charge_spin',help='Only Valid for `com\' file, default is `0 1\'',nargs='+')
    parser.add_argument('--gcheck',help='Only Valid for `com\' file, Gaussian check file name, default is input `com\' \
                        file name but with the extension `chk\'')
    parser.add_argument('--gtitle',help='Only Valid for `com\' file, default is `QM process\'',nargs='+')

    args = parser.parse_args()
    if args.basis_set is not None:
        line = ''
        for i in args.basis_set:
            line = line + i + ' '
        args.basis_set = line[:-1]

    if args.charge_spin is not None:
        line = ''
        for i in args.charge_spin:
            line = line + i + ' '
        args.charge_spin = line[:-1]

    if args.gtitle is not None:
        line = ''
        for i in args.gtitle:
            line = line + i + ' '
        args.gtitle = line[:-1]
        
            
    QMFile(args.file,args.atomnm,fextend=args.fextend,fname=args.fname,basis_set=args.basis_set,
           gcheck=args.gcheck,charge_spin=args.charge_spin,gtitle=args.gtitle).file_print()



