import numpy as np
from scipy import stats

from constants import L2MEAN, LCOV, D2MEAN, DCOV
from constants import GAMMA_D_ACI, GAMMA_L_ACI, GAMMA_D_GB, GAMMA_L_GB
from createRV.funcs import gblstats,lognstats
from pyStRe import ProbData, AnalysisOpt, Gfunc, CompReliab, SysReliab

def loadrv(coldesign, rhoLD=1.0, code='gb'):
    if code.lower() == 'gb':
        gammaD = GAMMA_D_GB
        gammaL = GAMMA_L_GB
    elif code.lower() == 'aci':
        gammaD = GAMMA_D_ACI
        gammaL = GAMMA_L_ACI
    Ntotal, Mtotal, e = coldesign.colcapacitymodel()
    # load Nd
    Ndnom = Ntotal/(gammaD+gammaL*rhoLD)
    Ndmean = Ndnom*D2MEAN
    Ndstd = Ndmean*DCOV
    rvNd = stats.norm(loc=Ndmean, scale=Ndstd)
    # load Nl
    Nlnom = rhoLD*Ntotal/(gammaD+gammaL*rhoLD)
    Nlmean = Nlnom*L2MEAN
    Nlstd = Nlmean*LCOV
    [locL, scaleL] = gblstats(Nlmean, Nlstd)
    rvNl = stats.gumbel_r(loc=locL, scale=scaleL)
    return rvNd, rvNl, Ndnom, Nlnom

def msr_form_6rv(rvnames, rvs, corr, e0=None, **kwargs):
    """rvnames= ['deltaN', 'deltae', 'N', 'M', 'D', 'L']
       rvs= coresponding rvs
       corr: matrix of coefficients of correlation of random variables
    """
    probdata = ProbData(names=rvnames, rvs=rvs, corr=corr, nataf=False)
    analysisopt = AnalysisOpt(gradflag='DDM', recordu=False, recordx=False,
            flagsens=False, verbose=False)
    # limit state 1
    def gf1(x, param=None):
        g = x[0]*x[2] - (x[4]+x[5])
        return g
    def dgdq1(x, param=None):
        dg0 = x[2]
        dg1 = 0.
        dg2 = x[0]
        dg3 = 0.
        dg4 = -1.
        dg5 = -1.
        return [dg0, dg1, dg2, dg3, dg4, dg5]
    gfunc1 = Gfunc(gf1, dgdq1)
    formBeta1 = CompReliab(probdata, gfunc1, analysisopt)
    # limit state 2
    def gf2(x, param=None, e0=e0):
        g = x[0]*x[3]*x[1] - (x[4]+x[5])*e0
        return g
    def dgdq2(x, param=None, e0=e0):
        dg0 = x[2]*x[1]
        dg1 = x[0]*x[3]
        dg2 = 0.
        dg3 = x[0]*x[1]
        dg4 = -e0
        dg5 = -e0
        return [dg0, dg1, dg2, dg3, dg4, dg5]
    gfunc2 = Gfunc(gf2, dgdq2)
    formBeta2 = CompReliab(probdata, gfunc2, analysisopt)

    # system reliability
    try:
        res1 = formBeta1.form_result()
        pf1 = res1.pf1
    except np.linalg.LinAlgError:
        pf1 = 0.
    try:
        res2 = formBeta2.form_result()
        pf2 = res2.pf1
    except np.linalg.LinAlgError:
        pf2 = 0.
    try:
        sysBeta = SysReliab([formBeta1, formBeta2], [-2])
        sysformres = sysBeta.mvn_msr(sysBeta.syscorr)
        pfsys = sysformres.pf
    except np.linalg.LinAlgError:
        pfsys = 0.
    return pfsys, pf1, pf2

def sysmc_6rv(rvnames, rvs, corr, e0=None, nsim=int(1e5)):
    """rvnames= ['deltaN', 'deltae', 'N', 'M', 'D', 'L']
       rvs= coresponding rvs
       corr: matrix of coefficients of correlation of random variables
       system Monte Carlo to verify msr_form
    """
    # model error samples
    meandNE = np.array([rvs[0].stats('m')[()], rvs[1].stats('m')[()]])
    stddNE = np.sqrt(np.array([[rvs[0].stats('v')[()], 0.],
        [0., rvs[1].stats('v')[()]]]))
    corrdNE = corr[:2, :2]
    covdNE = np.matmul(np.matmul(stddNE, corrdNE), stddNE)
    rvdNE = stats.multivariate_normal(mean=meandNE, cov=covdNE)
    dNEsmps = rvdNE.rvs(size=nsim)
    dNsmps = dNEsmps[:,0]
    dEsmps = dNEsmps[:,1]
    # resistance samples
    meanNM = np.array([rvs[2].stats('m')[()], rvs[3].stats('m')[()]])
    stdNM = np.sqrt(np.array([[rvs[2].stats('v')[()], 0.],
        [0., rvs[3].stats('v')[()]]]))
    corrNM = corr[-2:, -2:]
    covNM = np.matmul(np.matmul(stdNM, corrNM), stdNM)
    rvNM = stats.multivariate_normal(mean=meanNM, cov=covNM)
    NMsmps = rvNM.rvs(size=nsim)
    Nsmps = NMsmps[:,0]
    Msmps = NMsmps[:,1]
    # load samples
    Dsmps = rvs[4].rvs(size=nsim)
    Lsmps = rvs[5].rvs(size=nsim)
    # limit state 1
    gf1smps = dNsmps*Nsmps - (Dsmps+Lsmps)
    pf1 = np.sum(gf1smps<=0,dtype=float)/nsim
    # limit state 2
    gf2smps = dNsmps*Msmps*dEsmps - (Dsmps+Lsmps)*e0
    pf2 = np.sum(gf2smps<=0,dtype=float)/nsim
    # system reliability
    pfsys = 1.0 - np.sum((gf1smps>0)&(gf2smps>0), dtype=float)/nsim
    return pfsys, pf1, pf2

def msr_form(rvnames, rvs, corr, e0=None, **kwargs):
    """rvnames= ['N', 'M', 'D', 'L']
       rvs= coresponding rvs
       corr: matrix of coefficients of correlation of random variables
    """
    probdata = ProbData(names=rvnames, rvs=rvs, corr=corr, nataf=False)
    analysisopt = AnalysisOpt(gradflag='DDM', recordu=False, recordx=False,
            flagsens=False, verbose=False)
    # limit state 1
    def gf1(x, param=None):
        g = x[0] - (x[2]+x[3])
        return g
    def dgdq1(x, param=None):
        dg0 = 1.
        dg1 = 0.
        dg2 = -1.
        dg3 = -1.
        return [dg0, dg1, dg2, dg3]
    gfunc1 = Gfunc(gf1, dgdq1)
    formBeta1 = CompReliab(probdata, gfunc1, analysisopt)
    # limit state 2
    def gf2(x, param=None, e0=e0):
        g = x[1] - (x[2]+x[3])*e0
        return g
    def dgdq2(x, param=None, e0=e0):
        dg0 = 0.
        dg1 = 1.
        dg2 = -e0
        dg3 = -e0
        return [dg0, dg1, dg2, dg3]
    gfunc2 = Gfunc(gf2, dgdq2)
    formBeta2 = CompReliab(probdata, gfunc2, analysisopt)

    # system reliability
    try:
        res1 = formBeta1.form_result()
        pf1 = res1.pf1
    except np.linalg.LinAlgError:
        pf1 = 0.
    try:
        res2 = formBeta2.form_result()
        pf2 = res2.pf1
    except np.linalg.LinAlgError:
        pf2 = 0.
    try:
        sysBeta = SysReliab([formBeta1, formBeta2], [-2])
        sysformres = sysBeta.mvn_msr(sysBeta.syscorr)
        pfsys = sysformres.pf
    except np.linalg.LinAlgError:
        pfsys = 0.
    return pfsys, pf1, pf2

def sysmc(rvnames, rvs, corr, e0=None, nsim=int(1e5)):
    """rvnames= ['N', 'M', 'D', 'L']
       rvs= coresponding rvs
       corr: matrix of coefficients of correlation of random variables
       system Monte Carlo to verify msr_form
    """
    # resistance samples
    meanNM = np.array([rvs[0].stats('m')[()], rvs[1].stats('m')[()]])
    stdNM = np.sqrt(np.array([[rvs[0].stats('v')[()], 0.],
        [0., rvs[1].stats('v')[()]]]))
    [logNmean, logNstd] = lognstats(meanNM[0], stdNM[0,0])
    [logMmean, logMstd] = lognstats(meanNM[1], stdNM[1,1])
    corrNM = corr[-2:, -2:]
    covNM = np.matmul(np.matmul(stdNM, corrNM), stdNM)
    tmp = np.exp(logNmean+logMmean+0.5*(logNstd**2+logMstd**2))
    rNM = np.log(covNM[0,1]/tmp+1)
    meanlogNM = np.array([logNmean, logMmean])
    covlogNM = np.array([[logNstd**2, rNM], [rNM, logMstd**2]])
    rvlogNM = stats.multivariate_normal(mean=meanlogNM, cov=covlogNM)
    logNMsmps = rvlogNM.rvs(size=nsim)
    NMsmps = np.exp(logNMsmps)
    Nsmps = NMsmps[:,0]
    Msmps = NMsmps[:,1]
    # load samples
    Dsmps = rvs[2].rvs(size=nsim)
    Lsmps = rvs[3].rvs(size=nsim)
    # limit state 1
    gf1smps = Nsmps - (Dsmps+Lsmps)
    pf1 = np.sum(gf1smps<=0,dtype=float)/nsim
    # limit state 2
    gf2smps = Msmps - (Dsmps+Lsmps)*e0
    pf2 = np.sum(gf2smps<=0,dtype=float)/nsim
    # system reliability
    pfsys = 1.0 - np.sum((gf1smps>0)&(gf2smps>0), dtype=float)/nsim
    return pfsys, pf1, pf2

if __name__ == '__main__':
    import time
    import datetime
    from designcase import benchmark, loadcol
    from constants import DNMEAN, DNSTD, DEMEAN, DESTD, rhodNE, rhoNM
    from createRV.funcs import lognstats
    # resistance
    resistdata = np.load('./data/resistsmp.npz')
    resistmtr = resistdata['resistmtr']
    dNsmp = resistdata['dNsmp']
    dEsmp = resistdata['dEsmp']
    # N and M samples
    Nlhs = resistmtr[0,:]*dNsmp
    Mlhs = resistmtr[1,:]*dNsmp*dEsmp
    rhoNM4 = np.corrcoef(Nlhs, Mlhs)[0,1]
    Nmu4 = np.mean(Nlhs)
    Nsigma4 = np.std(Nlhs)
    Mmu4 = np.mean(Mlhs)
    Msigma4 = np.std(Mlhs)
    # rhoNM4 = 0.52548305
    # Nmu4 = 4458998.6649978179
    # Nsigma4 = 494197.8228842903
    # Mmu4 = 331677031.15943372
    # Msigma4 = 42903905.180368893

    np.random.seed(1)
    # for verification, model errors and resistances are assumed to be normal
    tfsol = np.load('./data/tfsol.npy')
    col0 = benchmark(tfsol)
    coldesign = loadcol(col0,code='gb', gammacc=0.9)
    [rvNd, rvNl, Ndnom, Nlnom] = loadrv(coldesign)
    # model errors
    [logDNmean, logDNstd] = lognstats(DNMEAN, DNSTD)
    # rvdN = stats.lognorm(logDNstd, scale=np.exp(logDNmean))
    rvdN = stats.norm(loc=DNMEAN, scale=DNSTD)
    [logDEmean, logDEstd] = lognstats(DEMEAN, DESTD)
    # rvdE = stats.lognorm(logDEstd, scale=np.exp(logDEmean))
    rvdE = stats.norm(loc=DEMEAN, scale=DESTD)
    # resistance smps
    resistdata = np.load('./data/resistsmp.npz')
    resistmtr = resistdata['resistmtr']
    Nmu = np.mean(resistmtr[0,:])
    Nsigma = np.std(resistmtr[0,:])
    [logNmean, logNstd] = lognstats(Nmu, Nsigma)
    # rvN = stats.lognorm(logNstd, scale=np.exp(logNmean))
    rvN = stats.norm(loc=Nmu, scale=Nsigma)
    Mmu = np.mean(resistmtr[1,:])
    Msigma = np.std(resistmtr[1,:])
    [logMmean, logMstd] = lognstats(Mmu, Msigma)
    # rvM = stats.lognorm(logNstd, scale=np.exp(logMmean))
    rvM = stats.norm(loc=Mmu, scale=Msigma)

    ## FORM + MSR (6 rvs)
    # reliability evaluation
    rvnames = ['deltaN', 'deltae', 'N', 'M', 'D', 'L']
    rvs = [rvdN, rvdE, rvN, rvM, rvNd, rvNl]
    corr = np.array([[1., rhodNE, 0., 0., 0., 0.],
                     [rhodNE, 1., 0., 0., 0., 0.],
                     [0., 0., 1., rhoNM, 0., 0.],
                     [0., 0., rhoNM, 1., 0., 0.],
                     [0., 0., 0., 0., 1., 0.],
                     [0., 0., 0., 0., 0., 1.]])
    start_delta_time = time.time()
    pfsys,pf1,pf2 = msr_form_6rv(rvnames, rvs, corr, col0.e0)
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    ## MC for verification (6rvs)
    start_delta_time = time.time()
    pfsysMC,pf1MC,pf2MC = sysmc_6rv(rvnames, rvs, corr, col0.e0, nsim=int(10e6))
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    ## FORM + MSR (4 rvs)
    # reliability evaluation
    # [logN4mean, logN4std] = lognstats(Nmu4, Nsigma4)
    # rvN4 = stats.lognorm(logN4std, scale=np.exp(logN4mean))
    # [logM4mean, logM4std] = lognstats(Mmu4, Msigma4)
    # rvM4 = stats.lognorm(logN4std, scale=np.exp(logM4mean))
    rvN4 = stats.norm(loc=Nmu4, scale=Nsigma4)
    rvM4 = stats.norm(loc=Mmu4, scale=Msigma4)
    rvnames4 = ['N', 'M', 'D', 'L']
    rvs4 = [rvN4, rvM4, rvNd, rvNl]
    corr4 = np.array([[1., rhoNM4, 0., 0.],
                     [rhoNM4, 1., 0., 0.],
                     [0., 0., 1., 0.],
                     [0., 0., 0., 1]])
    # e0 is magnified to increase pf
    start_delta_time = time.time()
    pfsys4,pf14,pf24 = msr_form(rvnames4, rvs4, corr4, col0.e0)
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    ## MC for verification (4rvs)
    start_delta_time = time.time()
    pfsysMC4,pf1MC4,pf2MC4 = sysmc(rvnames4, rvs4, corr4, col0.e0, nsim=int(10e6))
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
