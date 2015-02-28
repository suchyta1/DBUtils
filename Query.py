#!/usr/bin/env python

import desdb
import numpy as np
import esutil
import pyfits
import sys
import healpy as hp

from mpi4py import MPI
import mpifunctions
import DBfunctions

import numpy as np
import scipy.spatial
from sklearn.neighbors import NearestNeighbors

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



if __name__=='__main__': 

    select = {'table': 'sva1v2',
              'bands': ['i'],
              'truth': ['balrog_index', 'mag', 'ra', 'dec'],
              'sim': ['mag_auto']
             }


    cur = desdb.connect()
    if MPI.COMM_WORLD.Get_rank()==0:
        arr = cur.quick('SELECT unique(tilename) from balrog_%s_truth_%s' %(select['table'],select['bands'][0]), array=True)
        tiles = arr['tilename']
    else:
        tiles = None

    tiles = mpifunctions.Scatter(tiles)
    for tile in tiles:
        #truth = DBfunctions.TruthSelect(select, where="tilename='%s'"%(tile))
        #sim = DBfunctions.SimSelect(select, where="tilename='%s'"%(tile))
        #nosim = DBfunctions.NosimSelect(select, where="tilename='%s'"%(tile))
        
        where = "tilename='%s'"%(tile)
        truth, sim, nosim = DBfunctions.GetBalrog(select, where)
        print tile, len(truth), len(sim), len(nosim)
