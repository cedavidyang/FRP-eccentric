import numpy as np
import copy
from scipy import stats
import scipy.io as sio
from pyDOE import lhs

from designcase import benchmark, designcol, samplecol
from createRV.funcs import lognstats, wblstats

import time
import datetime
from pathos.multiprocessing import Pool


H2MEAN = 1.524; HSTD = 6.35
FCO2MEAN = 1.20; FCOCOV = 0.10
FS2MEAN = 1.20; FSCOV = 0.10
FFRP2MEAN = 1.25; FFRPCOV = 0.12
EFRP2MEAN = 1.0; EFRPCOV = 0.10
L2MEAN = 1.00; L2COV = 0.25
D2MEAN = 1.05; D2COV = 0.10

meNmean = 1.2176; meNcov = 0.10471
mee1mean = 1.85014; mee1cov = 0.17244
rNe1 = 0.6326

def samplecolarray(benchcol, nsim):
    # lhs samples
    lhdp = lhs(5, samples=nsim, criterion='corr', iterations=500)
    # random variable D
    hnom = benchcol.geo.h
    hmean = hnom + H2MEAN
    rvh = stats.norm(loc=hmean, scale=HSTD)
    hsmp = rvh.ppf(lhdp[:,0])
    # random variable fco
    fconom = benchcol.mat.fco
    fcomean = fconom*FCO2MEAN
    fcostd = fcomean*FCOCOV
    rvfco = stats.norm(loc=fconom, scale=fcostd)
    fcosmp = rvfco.ppf(lhdp[:,1])
    # random variable fs
    fsnom = benchcol.mat.fr
    fsmean = fsnom*FS2MEAN
    fscov = fsmean*FSCOV
    [logfsmean, logfsstd] = lognstats(fsmean, fscov)
    rvfs = stats.lognorm(logfsstd, scale=np.exp(logfsmean))
    fssmp = rvfs.ppf(lhdp[:,2])
    # random variable ffrp
    ffnom = benchcol.mat.fh
    ffmean = ffnom*FFRP2MEAN
    ffstd = ffmean*FFRPCOV
    [wblscale,wblc] = wblstats(ffmean, ffstd)
    rvff = stats.weibull_min(wblc, scale=wblscale)
    ffsmp = rvff.ppf(lhdp[:,3])
    # random variable Efrp
    Efrpnom = benchcol.mat.Eh
    Efrpmean = Efrpnom*EFRP2MEAN
    Efrpstd = Efrpmean*EFRPCOV
    [logEfmean, logEfstd] = lognstats(Efrpnom, Efrpstd)
    rvEf = stats.lognorm(logEfstd, scale=np.exp(logEfmean))
    Efsmp = rvEf.ppf(lhdp[:,4])
    # put into colArray
    colArray = np.empty((nsim,), dtype=object)
    for icol,col in enumerate(colArray):
        col = copy.deepcopy(benchcol)
        col = samplecol(col, h=hsmp[icol], fco=fcosmp[icol],
                fs=fssmp[icol], ff=ffsmp[icol], Ef=Efsmp[icol])
        colArray[icol] = col
    return colArray

def conductLHS(colArray,nprocess=1):
    start_delta_time = time.time()
    def mcevaluate(smp):
        return smp.colcapacitysection_fittedcurve()
    print 'CALC: Parallel version'
    try:
        pool = Pool(processes=nprocess)
        res = pool.map_async(mcevaluate, colArray).get(0xFFFF)
        pool.close()
        pool.join()
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
    return res

if __name__ == '__main__':
    np.random.seed(seed=1)
    col0 = benchmark()
    # N0,M0,et0 = col0.colcapacitysection_fittedcurve()
    nsim = 5
    colarraytest = samplecolarray(col0, nsim)
    # # series run
    # resistmtr = np.empty((2,nsim))
    # for i,col in enumerate(colarraytest):
        # Ni,Mi,ei = col.colcapacitysection_fittedcurve()
        # resistmtr[0,i] = Ni
        # resistmtr[1,i] = Mi
        # print 'MC: {} out of {}'.format(i+1, nsim)
    # parallel run
    res = conductLHS(colarraytest, nprocess=2)
    resistmtr = np.empty((3,nsim))
    for i,resi in enumerate(res):
        resistmtr[:,i] = res[i]
        print 'MC: {} out of {}'.format(i+1, nsim)

    # save results
    np.save('resistsmp.npy', resistmtr)
    sio.savemat('resistsmp.mat', {'resistmtr': resistmtr})
