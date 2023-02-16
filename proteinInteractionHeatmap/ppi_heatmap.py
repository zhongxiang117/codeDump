"""For protein-protein interaction analysis and the generation of heatmap

Declaration:
    Some of source codes are directly "copy-and-paste" from Schrodinger script.
    If you have license from them, then you are free to use those codes.
    Please use it at your own risks.
"""

import os
import sys

# hook to add needed modules
MMSHARE_EXEC = os.getenv('MMSHARE_EXEC')
if not MMSHARE_EXEC:
    raise ModuleNotFoundError('Codes should be executed as Schrodinger script')
sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(MMSHARE_EXEC)),'python/scripts')
)

from protein_interaction_analysis import getArgParser, validateOptions, run_job
from schrodinger import structure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#TODO:  axis label from "Residue Numbers" to "Residue ID"
FEATURES = [
    'version 0.1.0  : PPI-Heatmap, Feb 11th, 2023',
    'version 0.2.0  : make codes more compact',
    'version 0.3.0  : change `im` to `scatter`',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


if __name__ == '__main__':
    parser = getArgParser()
    opts = validateOptions(parser)
    run_job(opts)
    frame = pd.read_csv(opts.outfile)
    assert 'Residue' in frame.keys()
    assert 'Specific Interactions' in frame.keys()
    print(f'Note: results are saved to file: {opts.outfile}')

    chain_atomnum_offset = {}
    for ct in structure.StructureReader(opts.structure):
        tot = 0
        for c in ct.chain:
            chain_atomnum_offset[c.name] = tot
            tot += len([i for i in c.residue])
        chain_atomnum_offset['total'] = tot

    settings = {
        'default': {        # this key is important to indicate no interactions
            'key': 0,
            'color': 'white'
        },
        'hydrogen bond': {
            'key': 1,
            'color': 'red'
        },
        'salt bridge': {
            'key': 2,
            'color': 'yellow'
        },
        'pi stack': {
            'key': 3,
            'color': 'green'
        },
    }

    # important: replace all `nan` to empty string
    frame.fillna('', inplace=True)
    n = chain_atomnum_offset['total']
    defnum = settings['default']['key']
    table = [[defnum for i in range(n)] for j in range(n)]
    for i,j in zip(frame['Residue'],frame['Specific Interactions']):
        for v in j.split('\n'):
            bo = False
            for s in settings.keys():
                if s in v:
                    bo = True
                    vx = i.split(':')
                    x = chain_atomnum_offset[vx[0]] + int(vx[1])
                    vy = v.split()[-1].split(':')
                    y = chain_atomnum_offset[vy[0]] + int(vy[1])
                    table[x][y] = table[y][x] = settings[s]['key']
            if bo:
                break
    # reset self-interaction
    for i in range(n): table[i][i] = defnum

    table = np.array(table)
    row, col = np.where(table != defnum)
    ptable = table[sorted(list(set(row))),:]
    ptable = ptable[:,sorted(list(set(col)))]
    if len(ptable) == 0:
        print('Fatal: no interactions was found')
        sys.exit()

    fig, ax = plt.subplots()
    # construct (value,color) pairs to speed up calculation
    colors = {v['key']:v['color'] for v in settings.values()}
    settings.pop('default',None)
    for k in settings.keys():
        x,y = np.where(ptable == settings[k]['key'])
        ax.scatter(x,y,c=[colors[ptable[i][j]] for i,j in zip(x,y)],marker='s',label=k)
    ax.set_xlabel('Residue Numbers')
    ax.set_ylabel('Residue Numbers')
    ax.set_title('Protein Interaction Heatmap')
    ax.legend()
    plt.show()




