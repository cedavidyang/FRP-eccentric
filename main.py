import numpy as np
from scipy import stats
import scipy.io as sio
from scipy.optimize import minimize

from constants import *
from modelerror import testcolarray, test2model_ecc, test2model_Ru
from designcase import fcc2fco, benchmark, designcol, samplecol, loadcol
from resistanceLHS import samplecolarray, conductLHS
from reliability import msr_form, loadrv, sysmc
from createRV.funcs import gblstats, lognstats
from pyStRe import ProbData, AnalysisOpt, Gfunc, CompReliab, SysReliab

import time
import datetime
import os
import sys

if __name__ == '__main__':
    try:
        np.random.seed(2)
        nprocess = 4
        nlhs = 10000; iterations = 10000
        analysisNo = input('Reliability analysis number:')
        # some parameters
        if analysisNo == 1:
            frptype = 'cfrp'
            code = 'gb'
            colmodel = 'jiangteng13e1'
            fccmodel = 'tengetal09'
            eco = ECO_GB
            ecu = ECU_GB
        else:
            raise SystemExit('Analysis number not defined')

        # model error analysis
        print("MODEL ERROR ANALYSIS: BEGIN====================")
        if frptype.lower() in ['cfrp', 'c']:
            database = sio.loadmat('./database/frpdatabase.mat')['cfrpdatabase']
            if os.path.isfile('./data/modelerror-cfrp.mat'):
                medata = sio.loadmat('./data/modelerror-cfrp.mat')
                conductme = False
            else:
                conductme = True
        elif frptype.lower() in ['gfrp', 'g']:
            database = sio.loadmat('./database/frpdatabase.mat')['gfrpdatabase']
            if os.path.isfile('./data/modelerror-gfrp.mat'):
                medata = sio.loadmat('./data/modelerror-gfrp.mat')
                conductme = False
            else:
                conductme = True
        if conductme:
            colArray, Nu, Mu, e0 = testcolarray(database, failcode=1, fccmodel='tengetal09', code='gb')
            e1testArray, e1modelArray, e1sectionArray,\
                etestArray, emodelArray, esectionArray = test2model_ecc(colArray, Nu, Mu, e0, colmodel=colmodel, code='gb')
            NtestArray, NmodelArray, NsectionArray,\
                MtestArray, MmodelArray, MsectionArray = test2model_Ru(colArray, Nu, Mu, e0, colmodel=colmodel, code='gb')
            # save results
            sio.savemat('./data/modelerror-cfrp.mat', {'e1testArray': e1testArray,
                'e1modelArray': e1modelArray, 'e1sectionArray': e1sectionArray,
                'etestArray':etestArray, 'emodelArray':emodelArray,
                'esectionArray':esectionArray,
                'NsectionArray': NsectionArray, 'MsectionArray': MsectionArray,
                'NmodelArray': NmodelArray, 'MmodelArray': MmodelArray,
                'NtestArray': NtestArray, 'MtestArray': MtestArray
                })
        else:
            e1testArray = medata['e1testArray'].flatten()
            e1modelArray = medata['e1modelArray'].flatten()
            e1sectionArray = medata['e1sectionArray'].flatten()
            etestArray = medata['etestArray'].flatten()
            emodelArray = medata['emodelArray'].flatten()
            esectionArray = medata['esectionArray'].flatten()
            NtestArray = medata['NtestArray'].flatten()
            NmodelArray = medata['NmodelArray'].flatten()
            NsectionArray = medata['NsectionArray'].flatten()
            MtestArray = medata['MtestArray'].flatten()
            MmodelArray = medata['MmodelArray'].flatten()
            MsectionArray = medata['MsectionArray'].flatten()
        medict = {'dnmean': np.mean(NtestArray/NsectionArray),
                'dnstd': np.std(NtestArray/NsectionArray),
                'demean': np.mean(etestArray/esectionArray),
                'destd': np.std(etestArray/esectionArray)}
        medictmodel = {'dnmean': np.mean(NtestArray/NmodelArray),
                'dnstd': np.std(NtestArray/NmodelArray),
                'demean': np.mean(etestArray/emodelArray),
                'destd': np.std(etestArray/emodelArray)}
        print("MODEL ERROR ANALYSIS: END====================")


        print("REISISTANCE FITTING OF DESIGN CASES: BEGIN====================")
        # design cases (tf data)
        if os.path.isfile('./data/tfsol.npy'):
            tfsol = np.load('./data/tfsol.npy')
        else:
            nrho = FCC2FCOARRAY.shape[0]
            nfco = FCOARRAY.shape[0]
            tfsol = np.empty((nrho,nfco),dtype=float)
            for rhoindx, rhoi in enumerate(FCC2FCOARRAY):
                for fcoindx, fcoi in enumerate(FCOARRAY):
                    tfmin = 0.01*fcoi/eco*HDESIGN/(2*EHDESIGN)
                    sol = minimize(lambda x: abs(fcc2fco(x,fcoi)-rhoi), x0=1.0,
                            bounds=((tfmin,None),), options={'iprint':-1})
                    tf = sol.x[0]
                    tfsol[rhoindx, fcoindx] = tf
            np.save('./data/tfsol.npy', tfsol)
        # designed columns
        col0 = benchmark(tfsol)
        coldesignArray = [col0]
        colname = ['benchmark']
        for e2D in E2DARRAY:
            coli = designcol(tfsol, e2D=e2D)
            coldesignArray.append(coli)
            namei = 'e2D-'+str(e2D)
            colname.append(namei)
        for rhos in RHOSARRAY:
            coli = designcol(tfsol, rhos=rhos)
            coldesignArray.append(coli)
            namei = 'rhos-'+str(rhos)
            colname.append(namei)
        for d2D in DS2DARRAY:
            coli = designcol(tfsol, d2D=d2D)
            coldesignArray.append(coli)
            namei = 'd2D-'+str(d2D)
            colname.append(namei)
        for rhofcc in FCC2FCOARRAY:
            coli = designcol(tfsol, rhofcc=rhofcc)
            coldesignArray.append(coli)
            namei = 'rhofcc-'+str(rhofcc)
            colname.append(namei)
        for fco in FCOARRAY:
            coli = designcol(tfsol, fco=fco)
            coldesignArray.append(coli)
            namei = 'fco-'+str(fco)
            colname.append(namei)
        # generate resistance data of columns
        colfile = []
        lhsreslist = []
        for icol, coli in enumerate(coldesignArray):
            filenamei = './data/resistsmp-'+frptype+'-'+code+colname[icol]
            colfile.append(colfile)
            if os.path.isfile(filenamei+'.npz'):
                resistdata = np.load(filenamei+'.npz')
                lhsreslist.append(resistdata)
            else:
                print 'LHS: {} out of {} columns'.format(icol+1, np.shape(coldesignArray)[0])
                if 'e1' in colmodel:
                    colarraytest,dNsmp,dEsmp = samplecolarray(coli, nlhs,
                            medict=medictmodel, iterations=iterations)
                else:
                    colarraytest,dNsmp,dEsmp = samplecolarray(coli, nlhs,
                            medict=medict, iterations=iterations)
                # parallel run
                res = conductLHS(colarraytest, model=colmodel, nprocess=nprocess)
                resistmtr = np.empty((3,nlhs))
                for i,resi in enumerate(res):
                    resistmtr[:,i] = res[i]
                # save results
                np.savez(filenamei+'.npz', resistmtr=resistmtr, dNsmp=dNsmp, dEsmp=dEsmp)
                sio.savemat(filenamei+'.mat', {'resistmtr': resistmtr, 'dNsmp':dNsmp, 'dEsmp':dEsmp})
                # add to log
                lhsres = {'resistmtr': resistmtr, 'dNsmp':dNsmp, 'dEsmp':dEsmp}
                lhsreslist.append(lhsres)
        print("REISISTANCE FITTING OF DESIGN CASES: END====================")


        print("RELIABILITY-BASED CALIBRATION: BEGIN====================")
        # calibration based on reliability indices of different design cases
        pfsysArray = []
        pf1Array = []
        pf2Array = []
        betasysArray = []
        beta1Array = []
        beta2Array = []
        pfsysmcArray = []
        pf1mcArray = []
        pf2mcArray = []
        betasysmcArray = []
        beta1mcArray = []
        beta2mcArray = []
        gammaccArray = np.array([1.0])
        # for icol, coli in zip([16,17,18,19,20], np.array(coldesignArray)[16:21]):
        # gammaccArray = np.arange(0.75, 1.75, 0.05)
        for icol, coli in enumerate(coldesignArray):
            print 'CALIBRATION: {} out of {} columns'.format(icol+1, np.shape(coldesignArray)[0])
            # resistance
            resistdata = lhsreslist[icol]
            resistmtr = resistdata['resistmtr']
            dNsmp = resistdata['dNsmp']
            dEsmp = resistdata['dEsmp']
            # N and M samples
            Nlhs = resistmtr[0,:]*dNsmp
            Mlhs = resistmtr[1,:]*dNsmp*dEsmp
            rhoNM = np.corrcoef(Nlhs, Mlhs)[0,1]
            Nmu = np.mean(Nlhs)
            Nsigma = np.std(Nlhs)
            [logNmean, logNstd] = lognstats(Nmu, Nsigma)
            rvN = stats.lognorm(logNstd, scale=np.exp(logNmean))
            Mmu = np.mean(Mlhs)
            Msigma = np.std(Mlhs)
            [logMmean, logMstd] = lognstats(Mmu, Msigma)
            rvM = stats.lognorm(logMstd, scale=np.exp(logMmean))
            for gammacc in gammaccArray:
                # load effects
                coldesign = loadcol(coli,code=code, gammacc=gammacc)
                [rvNd, rvNl, Ndnom, Nlnom] = loadrv(coldesign, model='jiangteng13', code='gb')
                ## FORM + MSR (4 rvs)
                # reliability evaluation
                rvnames = ['N', 'M', 'D', 'L']
                rvs = [rvN, rvM, rvNd, rvNl]
                corr = np.array([[1., rhoNM, 0., 0.],
                                [rhoNM, 1., 0., 0.],
                                [0., 0., 1., 0.],
                                [0., 0., 0., 1]])
                # e0 is magnified to increase pf
                pfsys,pf1,pf2 = msr_form(rvnames, rvs, corr, coldesign.e0)
                beta1 = -stats.norm.ppf(pf1)
                beta2 = -stats.norm.ppf(pf2)
                betasys = -stats.norm.ppf(pfsys)
                pfsysmc,pf1mc,pf2mc = sysmc(rvnames, rvs, corr, coldesign.e0, nsim=int(5e6))
                beta1mc = -stats.norm.ppf(pf1mc)
                beta2mc = -stats.norm.ppf(pf2mc)
                betasysmc = -stats.norm.ppf(pfsysmc)
                # log results
                pfsysArray.append(pfsys)
                pf1Array.append(pf1)
                pf2Array.append(pf2)
                betasysArray.append(betasys)
                beta1Array.append(beta1)
                beta2Array.append(beta2)
                pfsysmcArray.append(pfsysmc)
                pf1mcArray.append(pf1mc)
                pf2mcArray.append(pf2mc)
                betasysmcArray.append(betasysmc)
                beta1mcArray.append(beta1mc)
                beta2mcArray.append(beta2mc)
        print("RELIABILITY-BASED CALIBRATION: END====================")


        # save data for postprocessing
        if frptype.lower() in ['cfrp', 'c']:
            filename = './data/calibration-cfrp.mat'
            filename_mc = './data/calibration-cfrp-mc.mat'
        elif frptype.lower() in ['gfrp', 'g']:
            filename = './data/calibration-gfrp.mat'
            filename_mc = './data/calibration-gfrp.mat'
        sio.savemat(filename, {'pfsysArray': pfsysArray, 'pf1Array':pf1Array, 'pf2Array':pf2Array,
            'betasysArray':betasysArray, 'beta1Array':beta1Array, 'beta2Array':beta2Array})
        sio.savemat(filename_mc, {'pfsysArray': pfsysmcArray, 'pf1Array':pf1mcArray, 'pf2Array':pf2mcArray,
            'betasysArray':betasysmcArray, 'beta1Array':beta1mcArray, 'beta2Array':beta2mcArray})
    finally:
        from emailreminder import send_notification
        # send_notification()
