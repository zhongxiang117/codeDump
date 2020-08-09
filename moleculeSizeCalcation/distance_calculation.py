#!/usr/bin/env python3

import argparse

class Size_calculation(object):
    def __init__(self,pdbfile):
        self.pdbfile = pdbfile
        self.totlist = self._pro_pdbfile()

    def run(self):
        if self.totlist == 0:
            print('Error: the PDB file is wrong, no inputs')
            exit()
        elif self.totlist == 1:
            print('Warning: the input molecule only contains one atom!')
            exit()
        else:
            rmax_1d,(i,j) = self.distance_1D()
            print('\nThe input molecule\'s 1D longest distance is: %.3f Angstrom' % rmax_1d)
            print('The corresponded two atoms\' indices are: (%d, %d)' % (i+1,j+1))

            rmax_2d,(p,w) = self.distance_2D()
            print('\nThe input molecule\'s 2D longest distance is: %.3f Angstrom' % rmax_2d)
            print('The corresponded two atoms\' indices are: (%d, %d)\n' % (p+1,w+1))


    def _pro_pdbfile(self):
        totlist = []
        with open(self.pdbfile,mode='rt') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break
                else:
                    xyz = []
                    if len(line) >= 6 and line[:3].upper() == 'END':
                            break
                    elif len(line) >= 54:
                        xyz.append(float(line[29:38]))
                        xyz.append(float(line[38:46]))
                        xyz.append(float(line[46:54]))
                        totlist.append(xyz)
        return totlist


    def distance_1D(self):
        '''first dimension of longest distance between two atoms'''
        # the maximum distance
        rmax = 0
        # the index of atoms
        indx = jndx = 0
        i = 1
        while i < len(self.totlist):
            j = 0
            while j < i:
                x = self.totlist[j][0] - self.totlist[i][0]
                y = self.totlist[j][1] - self.totlist[i][1]
                z = self.totlist[j][2] - self.totlist[i][2]
                rcalc = x*x + y*y + z*z
                if rcalc >= rmax:
                    rmax = rcalc
                    indx = i
                    jndx = j
                j += 1
            i += 1
        return pow(rmax,0.5), (indx,jndx)


    def distance_2D(self):
        '''second dimension of longest distance between two atoms,
           which is perpendicular with the first dimension distance's two atoms'''
        r,(indx,jndx) = self.distance_1D()

        s = self.totlist[indx]
        t = self.totlist[jndx]

        stx = s[0] - t[0]
        sty = s[1] - t[1]
        stz = s[2] - t[2]

        # make a copy of total coordinates to avoid any accidentally modifications
        ndxlist = self.totlist

        # the maximum distance
        rmax = 0
        # the index of atoms
        self.pndx = self.wndx = 0
        i = 1
        while i < len(ndxlist):
            j = 0
            while j < i:
                m = ndxlist[j][0]*stx + ndxlist[j][1]*sty + ndxlist[j][2]*stz
                ssum = s[0]*stx + s[1]*sty + s[2]*stz
                tsum = t[0]*stx + t[1]*sty + t[2]*stz

                delta_sm = m - ssum
                delta_tm = m - tsum

                # take care of the same atom as well as the round-off-error
                if round(delta_sm,6) == 0 or round(delta_tm,6) == 0:
                    rcalc = 0
                else:
                    u = delta_sm / delta_tm

                    # get the projected-point coordinates
                    x = (s[0] - t[0]*u) / (1-u)
                    y = (s[1] - t[1]*u) / (1-u)
                    z = (s[2] - t[2]*u) / (1-u)

                    # calculate the perpendicular vector and its modulus
                    vx = ndxlist[j][0] - x
                    vy = ndxlist[j][1] - y
                    vz = ndxlist[j][2] - z
                    v_modulus = pow(vx*vx + vy*vy + vz*vz, 0.5)

                    # calculate two-points base vector
                    bx = ndxlist[j][0] - ndxlist[i][0]
                    by = ndxlist[j][1] - ndxlist[i][1]
                    bz = ndxlist[j][2] - ndxlist[i][2]

                    # calculate the final 2D distance
                    rcalc = abs( (bx*vx + by*vy + bz*vz) / v_modulus )

                if rcalc >= rmax:
                    rmax = rcalc
                    self.pndx = i
                    self.wndx = j

                j += 1
            i += 1
        return rmax, (self.pndx,self.wndx)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='molecule_size_calculation')
    parser.add_argument('-f','--pdb',help='input pdb file',required=True)
    args = parser.parse_args()
    Size_calculation(args.pdb).run()


