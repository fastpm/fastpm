#!/usr/bin/python3

import numpy as np
import nbodykit.lab as nkit
from nbodykit.source.catalog import BigFileCatalog
import os

def compute_pk(base,Nmesh,BoxSize=250.):
    headerdir = base+'Header/'
    datadir   = base+'1/'

    kmin = 2*np.pi/BoxSize
    dk   = kmin
    kmax = Nmesh*np.pi/BoxSize

    cat  = BigFileCatalog(datadir,header=headerdir)
    mesh = cat.to_mesh(BoxSize=BoxSize,Nmesh=Nmesh,resampler='cic',compensated=True)
    results = nkit.FFTPower(mesh,mode='1d',kmin=kmin,dk=dk,kmax=kmax)
    P  = results.power
    k  = P['k']
    Pk = P['power'].real
    return k,Pk


base  = '../../output_tests/snaps/'
baseout  = 'output_pk/'
Nmesh = 2048
snaps = os.listdir(base)
for i,snap in enumerate(snaps):
    print(snap)
    if snap !='old':
        fin    = base+snap+'/'
        k,Pk   = compute_pk(fin,Nmesh)
        tosave = np.array([k,Pk]).T
        fout   = baseout+snap+'.dat'
        np.savetxt(fout,tosave)
    else:
        print('Hello there!')
    print('Done step ',i+1,' / ',len(snaps))
print('Done!')

