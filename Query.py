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
              'des': 'sva1_coadd_objects',
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
        
        where = "where tilename='%s'"%(tile)
        truth, sim, nosim = DBfunctions.GetBalrog(select, truthwhere=where, simwhere=where)
        des = DBfunctions.GetDES(select, where=where)
        print tile, len(truth), len(sim), len(nosim), len(des)

    """
    nside = 64
    balrog_nest = False
    if MPI.COMM_WORLD.Get_rank()==0:
        print hp.nside2npix(nside)
        dbrange = DBfunctions.GetCourseRange(select)
        ramin, ramax, decmin, decmax, hpindex = DBfunctions.GetHealPixRectangles(nside, dbrange, balrog_nest)
    else:
        ramin, ramax, decmin, decmax, hpindex = [None]*5

    ramin, ramax, decmin, decmax, hpindex = mpifunctions.Scatter(ramin, ramax, decmin, decmax, hpindex)
    nonzero = DBfunctions.DoCount(ramin, ramax, decmin, decmax, hpindex, select)
    ramin, ramax, decmin, decmax = mpifunctions.Gather(ramin, ramax, decmin, decmax)
    hpindex = mpifunctions.Gather(hpindex, dtype=np.int64)
    nonzero = mpifunctions.Gather(nonzero, dtype=np.bool_)

    if MPI.COMM_WORLD.Get_rank()==0:
        ramin, ramax, decmin, decmax, hpindex = DBfunctions.Cut(ramin, ramax, decmin, decmax, hpindex, cut=nonzero)
        print hpindex.shape, hp.nside2pixarea(nside, degrees=True)

    ramin, ramax, decmin, decmax, hpindex = mpifunctions.Scatter(ramin, ramax, decmin, decmax, hpindex)
    for i in range(len(hpindex)):
        twhere = 'WHERE dec between %f and %f and ra between %f and %f' %(decmin[i],decmax[i], ramin[i],ramax[i])
        swhere = 'WHERE deltawin_j2000 between %f and %f and alphawin_j2000 between %f and %f' %(decmin[i],decmax[i], ramin[i],ramax[i])
        dws = []
        for band in select['bands']:
            dws.append('deltawin_j2000_%s between %f and %f and alphawin_j2000_%s between %f and %f' %(band,decmin[i],decmax[i], band,ramin[i],ramax[i]) )
        dwhere = 'WHERE %s' %(' and '.join(dws))

        truth, sim, nosim = DBfunctions.GetBalrog(select, truthwhere='', simwhere=swhere)
        des = DBfunctions.GetDES(select, where=dwhere)
        print hpindex[i], len(truth), len(sim), len(nosim)
    """
