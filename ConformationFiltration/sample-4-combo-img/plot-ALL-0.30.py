#!/usr/bin/env python3

import os
import sys
import argparse
import matplotlib.pyplot as plt


FEATURES = [
    'version 0.10 : start',
    'version 0.20 : split to functions',
    'version 0.30 : add argparse RELEASE',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


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




def readfile(file):
    """read file: bulk-probability-data.txt

    Return:

        filelist    :   1D List     :   [file, file, ...]

        plotdatalist:   dict
            bonds:  List[ List[List[float]] ], }
                => None, means not exit
                => [ [ini-X, ini-Y, fin-X, fin-Y ] ]
            angles:  same as bonds

        tolerance   :  1D List      :   [bonds, angles]
    
    Note:
        special file < OVERALL > is ignored
    """
    with open(file,'rt') as f:
        profile = f.readlines()

    i = 0
    tolerance = None
    filelist = []
    datalist = []
    while i < len(profile):
        line = profile[i]
        line = line[:line.find('#')] if line.find('#') != -1 else line
        line = line.strip()
        if len(line) == 0:
            i += 1
            continue

        if line.find('@TOLERANCE') != -1:
            tmp = line.split()
            tolerance = [float(tmp[1]), float(tmp[2])]
            i += 1
            continue
        
        elif line.find('@FILE') != -1:
            filelist.append(line.split()[1])

            j = i + 1
            fdict = {'bonds':[], 'angles':[]}
            while j < len(profile):
                line = profile[j]
                line = line[:line.find('#')] if line.find('#') != -1 else line
                line = line.strip()
                if line.find('@FILE') != -1:
                    break
                if len(line) == 0:
                    j += 1
                    continue
                if line.find('@BONDS') != -1 and line.find('ALL') != -1:
                    fdict['bonds'].append(line)
                elif line.find('@ANGLES') != -1 and line.find('ALL') != -1:
                    fdict['angles'].append(line)
                j += 1
            datalist.append(fdict)
            i = j
        else:
            i += 1

    if tolerance is None: tolerance = [0.1, 0.1]

    for cnt,t in enumerate(datalist):
        if len(t['bonds']) % 4 != 0 or len(t['angles']) % 4 != 0:
            print('Warning: wrong input file < {:} >'.format(file))
            print('  => for processing file < {:} >'.format(filelist[cnt]))
            raise ValueError('wrong input')

    # format  :  dict
    #   key:
    #       bonds:  List[ List[List[float]] ], }
    #           => None, means not exit
    #           => [ [ini-X, ini-Y, fin-X, fin-Y ] ]
    #       angles
    plotdatalist = {'bonds':[], 'angles':[]}
    dtb = tolerance[0]
    dta = tolerance[1]
    for data in datalist:
        i = 3
        while i < len(data['bonds']):
            inix = data['bonds'][i-3]
            iniy = data['bonds'][i-2]
            yil = iniy.split()
            if len(yil) <= 5:
                i += 4
                plotdatalist['bonds'].append([None,None,None,None])
                continue
            yil = [int(s) for s in yil[1:]]
            v = float(inix.split()[1])
            xil = [v+dtb*s for s in range(len(yil))]

            finx = data['bonds'][i-1]
            finy = data['bonds'][i]
            yfl = finy.split()
            if len(yfl) <= 5:
                i += 4
                plotdatalist['bonds'].append([None,None,None,None])
                continue
            yfl = [int(s) for s in yfl[1:]]
            v = float(finx.split()[1])
            xfl = [v+dtb*s for s in range(len(yfl))]
            plotdatalist['bonds'].append([xil,yil,xfl,yfl])
            i += 4
        
        i = 3
        while i < len(data['angles']):
            inix = data['angles'][i-3]
            iniy = data['angles'][i-2]
            yil = iniy.split()
            if len(yil) <= 5:
                i += 4
                plotdatalist['angles'].append([None,None,None,None])
                continue
            yil = [int(s) for s in yil[1:]]
            v = float(inix.split()[1])
            xil = [v+dta*s for s in range(len(yil))]

            finx = data['angles'][i-1]
            finy = data['angles'][i]
            yfl = finy.split()
            if len(yfl) <= 5:
                i += 4
                plotdatalist['angles'].append([None,None,None,None])
                continue
            yfl = [int(s) for s in yfl[1:]]
            v = float(finx.split()[1])
            xfl = [v+dta*s for s in range(len(yfl))]
            plotdatalist['angles'].append([xil,yil,xfl,yfl])
            i += 4

    return filelist,plotdatalist,tolerance




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




def save_image(datalist,key='bonds',dt=0.1,fname=None):
    # the number of linestyle should be in ODD number
    linestyle = CYCLE(['solid','dotted','dashdot',])
    colors = CYCLE(['b','r','m','g','y','brown','palegreen','deepskyblue'])
    if key.lower() == 'bonds':
        info = 'Filtration on Bonds ({:} Angstrom)'.format(dt)
    else:
        info = 'Filtration on Angles ({:} Degree)'.format(dt)

    if fname is None: fname = 'bulk-image-compare-{:}.jpg'.format(key)

    # format 2D : List[ List[int,int] ]
    datalist = [i for i in datalist if i[0] is not None]
    
    molnms = [[sum(i[1]),sum(i[3])] for i in datalist]
    reflist = [i[0] for i in molnms]
    reflist = sorted(range(len(reflist)), key=lambda k: reflist[k])

    sortdatalist = [datalist[i] for i in reflist]

    lines = []
    for d in sortdatalist:
        color = colors.next()
        ils = linestyle.next()
        fls = linestyle.next()
        itxt = 'Initial Mols = {:}'.format(sum(d[1]))
        ftxt = 'Final Mols = {:}'.format(sum(d[3]))
        iln, = plt.plot(d[0],d[1],color=color,linestyle=ils,label=itxt)
        fln, = plt.plot(d[2],d[3],color=color,linestyle=fls,label=ftxt)
        lines.append(iln)
        lines.append(fln)
    plt.legend(handles=lines,loc='upper right')
    plt.title(info)
    print('Note: image file is saved to < {:} >'.format(fname))
    plt.savefig(fname)
    plt.close()




def parsecmd():
    """Parse command line input"""
    parser = argparse.ArgumentParser(
        description='Plot compare results from Confromation Filtration',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION,
    )
    parser.add_argument(
        '-f','--file',
        help='File: bulk-probability-data',
    )
    parser.add_argument(
        '--features',
        help='Show develop features',
        action='store_true',
    )
    parser.add_argument(
        '-o','--fname',
        help='Output file name',
    )

    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    args = parser.parse_args(sys.argv[1:])
    if args.features:
        for i in FEATURES:
            print(i)
        exit()
    
    if args.file is None:
        print('Warning: -f/--file is missing')
        exit()

    return args




def main():
    args = parsecmd()

    filelist,plotdatalist,tolerance = readfile(args.file)

    if len(filelist) == 0:
        print('Warning: No inputs')
        exit()

    boall = False if filelist[-1].find('OVERALL') == -1 else True

    bondlist = plotdatalist['bonds']
    if boall: bondlist = bondlist[:-1]
    if args.fname is None:
        fname = 'bulk-image-compare-bonds' 
        fname = file_gen_new(fname,'jpg')
    else:
        fname = file_gen_new(args.fname,'jpg')
    save_image(bondlist,key='bonds',dt=tolerance[0],fname=fname)


    anglelist = plotdatalist['angles']
    if boall: anglelist = anglelist[:-1]
    if args.fname is None:
        fname = 'bulk-image-compare-angles' 
        fname = file_gen_new(fname,'jpg')
    else:
        fname = file_gen_new(args.fname,'jpg')
    save_image(anglelist,key='angles',dt=tolerance[1],fname=fname)




if __name__ == '__main__':
    main()