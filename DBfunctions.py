import os
import sys
import subprocess
import copy
import copy
import desdb
import numpy as np
import esutil
import pyfits
from mpi4py import MPI


def AngularDistance(r1, r2, d1, d2):
    ra1, ra2, dec1, dec2 = deg2rad(r1, r2, d1, d2)
    ddec = np.radians(dec2-dec1)
    dra = np.radians(ra2-ra1)
    a = np.sin(ddec/2) * np.sin(ddec/2) + np.cos(dec1) * np.cos(dec2) * np.sin(dra/2) * np.sin(dra/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return np.degrees(c)


def QuerySFW(s, f, w, select, where):
    cur = desdb.connect()
    if where != None:
        for i in range(len(select['bands'])):
            ww = 'truth_%s.%s' %(select['bands'][i],where)
            w.insert(0, ww)

    selects = 'SELECT %s' %(', '.join(s))
    froms = 'FROM %s' %(', '.join(f))
    wheres = 'WHERE %s' %(' and '.join(w))
    q = '%s %s %s' %(selects, froms, wheres)
    return cur.quick(q, array=True)


def TruthSelect(select, where=None):
    f = []
    w = []
    s = []
    band0 = select['bands'][0]
    for i in range(len(select['bands'])):
        band = select['bands'][i]
        f.append('balrog_%s_truth_%s truth_%s' %(select['table'], band, band) )
        for sel in select['truth']:
            s.append('truth_%s.%s as %s_%s' %(band,sel,sel,band))
        if i > 0:
            w.append('truth_%s.balrog_index=truth_%s.balrog_index' %(band0,band) )

    return QuerySFW(s, f, w, select, where)


def SimSelect(select, where=None):
    f = []
    w = []
    s = []
    band0 = select['bands'][0]

    '''
    qq = GetQuery(select, where, kind='truth')
    cur = desdb.connect()
    a = cur.quick(qq, array=True)
    print len(a), a.dtype.names
    sys.exit()
    '''

    for i in range(len(select['bands'])):
        band = select['bands'][i]

        f.append('balrog_%s_sim_%s sim_%s' %(select['table'], band, band) )
        f.append('balrog_%s_truth_%s truth_%s' %(select['table'], band, band) )
        for sel in select['sim']:
            s.append('sim_%s.%s as %s_%s' %(band,sel,sel,band))
        for sel in select['truth']:
            s.append('truth_%s.%s as %s_%s' %(band,sel,sel,band))
        w.append('truth_%s.balrog_index=sim_%s.balrog_index' %(band,band) )
        if i > 0:
            w.append('truth_%s.balrog_index=truth_%s.balrog_index' %(band0,band) )

    return QuerySFW(s, f, w, select, where)


def NosimSelect(select, where=None):
    f = []
    w = []
    s = []

    band0 = select['bands'][0]
    f.append('balrog_%s_nosim_det nosim_det' %(select['table']) )
    for sel in select['sim']:
        s.append('nosim_det.%s as %s_det' %(sel,sel))
    for i in range(len(select['bands'])):
        band = select['bands'][i]
        f.append('balrog_%s_truth_%s truth_%s' %(select['table'], band, band) )
        for sel in select['truth']:
            s.append('truth_%s.%s as %s_%s' %(band,sel,sel,band))
        w.append('truth_%s.balrog_index=nosim_det.balrog_index' %(band) )

    return QuerySFW(s, f, w, select, where)


def AppendSelects(ss, select, truthband=None, simband=None, truthlabel='truth', simlabel='sim', simbandappend='', truthbandappend=''):
    for sel in select['truth']:
        ss.append('%s_%s.%s as %s%s' %(truthlabel,truthband,sel,sel,truthbandappend))

    if simlabel != 'truth':
        for sel in select['sim']:
            ss.append('%s_%s.%s as %s%s' %(simlabel,simband,sel,sel,simbandappend))
    return ss


def GetQuery(select, where, kind='sim'):
    js = []
    ss = []
    bands = select['bands']
    for i in range(len(bands)):
        truthband = bands[i]
        simband = bands[i]
        truthbandappend = '_%s' %(bands[i])
        simbandappend = '_%s' %(bands[i])
        simlabel = 'joined'
        if kind=='nosim':
            simband = 'det'
            if i==0:
                simbandappend = '_det'
            else:
                simlabel = 'truth'
        elif kind=='truth':
            simlabel = 'truth'

        bq = Join2Truth(select, bands[i], simband, where, kind=kind)
        ss = AppendSelects(ss, select, truthband=bands[i], simband=bands[i], truthlabel='joined', simlabel=simlabel, truthbandappend=truthbandappend, simbandappend=simbandappend)
        if i==0:
            f = """FROM (%s) joined_%s"""%(bq, bands[i])
        else:
            js.append("""JOIN (%s) joined_%s ON joined_%s.balrog_index=joined_%s.balrog_index"""%(bq, bands[i], bands[i], bands[0]))

    q = 'SELECT %s %s %s' %(', '.join(ss), f, ' '.join(js))
    return q


def Join2Truth(select, truthband, simband, where, kind='sim'):
    ss = []
    ss = AppendSelects(ss, select, truthband=truthband, simband=simband, truthlabel='truth', simlabel=kind)
    ss = ', '.join(ss)
    f = """(SELECT %s FROM balrog_%s_truth_%s where %s) truth_%s"""%(', '.join(select['truth']), select['table'],truthband, where, truthband)
    q = """SELECT %s FROM %s""" %(ss, f)
    if kind in ['sim','nosim']:
        j = """(SELECT %s, balrog_index FROM balrog_%s_%s_%s where %s) %s_%s"""%(', '.join(select['sim']), select['table'], kind, simband, where, kind, simband)
        o = """truth_%s.balrog_index=%s_%s.balrog_index""" %(truthband,kind,simband)
        q = """%s JOIN %s ON %s"""%(q,j,o)
    return q


def GetBalrog(select, where):
    cur = desdb.connect()

    q = GetQuery(select, where, kind='truth')
    truth = cur.quick(q, array=True)

    q = GetQuery(select, where, kind='sim')
    sim = cur.quick(q, array=True)

    q = GetQuery(select, where, kind='nosim')
    nosim = cur.quick(q, array=True)
    
    return [truth, sim, nosim]