import sys
import numpy as np
import scipy.io as sio
from scipy.optimize import minimize
from scipy.interpolate import interp1d

from section import Material, Geometry
from column import Column
from constants import *

def testcolarray(database, failcode=1, **kwargs):
    """kwargs should at least include fccmodel and code"""
    colArray = []
    failmode = np.array([columnData[FMCOL] for columnData in database])
    NuArray = np.array([columnData[NUCOL]*1e3 for columnData in database])
    MuArray = np.array([columnData[MUCOL]*1e6 for columnData in database])
    e0Array = np.array([columnData[E0COL] for columnData in database])
    NuArray = NuArray[failmode==failcode]
    MuArray = MuArray[failmode==failcode]
    e0Array = e0Array[failmode==failcode]
    # concrete parameters
    code = kwargs['code']
    if code.lower() == 'gb':
        ecu = ECU_GB
        eco = ECO_GB
    else:
        raise SystemExit('Design code not available (for now)')
    for icol, columnData in enumerate(database[failmode==failcode]):
        # define material
        fc = columnData[FCCOL]
        fr = columnData[FRCOL]
        Er = columnData[ERCOL]*1e3
        # rupture strain of steel is assmued to be 100 times larger than yield strain
        eru = fr/Er*100
        # reinss must be a vector function
        reinss = lambda x,Er=Er,fr=fr: np.maximum(np.minimum(x*Er, fr), -fr)
        mat = Material(fc, fr, eco, ecu, eru, Er=Er, concss=None, reinss=reinss)
        eh = columnData[EHCOUPON]
        Eh = columnData[EHCOL]*1e3
        fh = eh*Eh
        mat.setconfinemat(fh, Eh, eh)
        # define geometry
        h = columnData[HCOL]
        Ac = np.pi*h**2/4.
        Afb = np.array([0.])    #frp bars
        tft = columnData[TFTCOL]
        Astotal = columnData[RHOSCOL]*Ac
        ns = columnData[NLCOL]
        tmp = np.arange(1., ns+1., 2)
        nsportion = np.insert(tmp, [0, ns/2], [0., ns])
        Asportion = (nsportion[1:]-nsportion[:-1])/ns
        As = Astotal*Asportion
        xf = np.array([0.])    #frp bars
        # steel location
        alphas = np.linspace(0., 1., ns/2+1)*np.pi
        rs = h/2. - columnData[COVERCOL]
        xs = h/2. - rs*np.cos(alphas)
        def bdist(x, h=h, b=h):
            if x<0 or x>h:
                bsec = 0.
            else:
                bsec = 2*np.sqrt((h/2)**2-(x-h/2)**2)
            return bsec
        geo = Geometry('rc', h, Ac, Afb, tft, As, xf, xs, bdist)
        # define column
        column = Column(geo=geo, mat=mat, cost=cost)
        column.setfrpcolmat(model=kwargs['fccmodel'])
        l = columnData[L2HCOL]*columnData[HCOL]
        e0 = e0Array[icol]
        column.setcolgeo(l, e0)
        colArray.append(column)
    return colArray, NuArray, MuArray, e0Array

def test2model_ecc(colArray, NuArray, MuArray, e0Array, method=None, **kwargs):
    """kwargs should at least include colmodel"""
    # extra eccentricity
    e1modelArray = []
    e1sectionArray = []
    e1testArray = []
    # total eccentricity
    emodelArray = []
    esectionArray = []
    etestArray = []
    for icol,column in enumerate(colArray):
        Mu = MuArray[icol]
        Nu=NuArray[icol]
        e0 = e0Array[icol]
        etotal = Mu/Nu
        e1test = etotal - e0
        if method is None:
            Ncolmodel, Mcolmodel, e1model = column.colcapacitymodel(model=kwargs['colmodel'])
            print 'model: {} out of {}'.format(icol+1, np.shape(colArray)[0])
            Ncolsection, Mcolsection, e1section = column.colcapacitysection_fittedcurve()
            print 'section: {} out of {}'.format(icol+1, np.shape(colArray)[0])
            # add to log
            e1testArray.append(e1test)
            e1modelArray.append(e1model)
            e1sectionArray.append(e1section)
            etestArray.append(etotal)
            emodelArray.append(Mcolmodel/Ncolmodel)
            esectionArray.append(Mcolsection/Ncolsection)
        elif method.lower() == 'section':
            Ncolsection, Mcolsection, e1section = column.colcapacitysection_fittedcurve()
            print 'section: {} out of {}'.format(icol+1, np.shape(database[failmode==1])[0])
            # add to log
            e1testArray.append(e1test)
            e1sectionArray.append(e1section)
            etestArray.append(etotal)
            esectionArray.append(Mcolsection/Ncolsection)
        elif method.lower() == 'model':
            Ncolmodel, Mcolmodel, e1model = column.colcapacitymodel(model=kwargs['colmodel'])
            print 'model: {} out of {}'.format(icol+1, np.shape(database[failmode==1])[0])
            # add to log
            e1testArray.append(e1test)
            e1modelArray.append(e1model)
            etestArray.append(etotal)
            emodelArray.append(Mcolmodel/Ncolmodel)
    e1testArray = np.array(e1testArray)
    e1modelArray = np.array(e1modelArray)
    e1sectionArray = np.array(e1sectionArray)
    etestArray = np.array(etestArray)
    emodelArray = np.array(emodelArray)
    esectionArray = np.array(esectionArray)
    return e1testArray, e1modelArray, e1sectionArray, etestArray, emodelArray, esectionArray

def test2model_Ru(colArray, NuArray, MuArray, e0Array, method=None, saveres=True, **kwargs):
    """kwargs should at least include code"""
    # ultimate capacity
    MsectionArray = []
    NsectionArray = []
    MmodelArray = []
    NmodelArray = []
    MtestArray = []
    NtestArray = []
    for icol,column in enumerate(colArray):
        Mu = MuArray[icol]
        Nu=NuArray[icol]
        e0 = e0Array[icol]
        etotal = Mu/Nu
        e1test = etotal - e0
        if method is None or method.lower() == 'section':
            # section analysis
            Nccu = column.axial_compression(confine='yes')
            NarraySection = np.linspace(1e-3, 1-1e-3, 50)*Nccu
            MarraySection = []
            ciarray = []
            for Ni in NarraySection:
                Nsec,Msec,ci = column.capacity(Ni)
                MarraySection.append(Msec)
                ciarray.append(ci)
            MarraySection = np.array(MarraySection)
            fNM = interp1d(NarraySection, MarraySection, bounds_error=False,
                    fill_value=(MarraySection[0], MarraySection[-1]))
            def desection(Niscale, scale, fNM=fNM, etotal=etotal):
                Ni = Niscale*scale
                Msec = fNM(Ni)[()]
                return abs(Msec/Ni-etotal)
            scale = 1e3
            sol1 = minimize(desection, Nu/scale, args=(scale,),
                    bounds=((0,Nccu/scale),), options={'iprint':-1})
            Nsection = sol1.x[0]*scale
            Msection = Nsection*etotal
            print 'section: {} out of {}'.format(icol+1, np.shape(colArray)[0])
        if method is None or method.lower() == 'model':
            # jiang and teng (2013) model
            def demodel(thetai, column=column, etotal=etotal):
                Nmodel = column.modelNmat(thetai)
                Mmodel = column.modelMmat(thetai)
                return abs(Mmodel/Nmodel-etotal)
            sol2 = minimize(demodel, (0.1+0.75)/2., bounds=((0.1,0.75),))
            theta = sol2.x[0]
            Nmodel = column.modelNmat(theta)
            Mmodel = column.modelMmat(theta)
            print 'model: {} out of {}'.format(icol+1, np.shape(colArray)[0])
        # add to log
        if method is None:
            NtestArray.append(Nu)
            MtestArray.append(Mu)
            NsectionArray.append(Nsection)
            MsectionArray.append(Msection)
            NmodelArray.append(Nmodel)
            MmodelArray.append(Mmodel)
        elif method.lower() == 'section':
            NtestArray.append(Nu)
            MtestArray.append(Mu)
            NsectionArray.append(Nsection)
            MsectionArray.append(Msection)
        elif method.lower() == 'model':
            NtestArray.append(Nu)
            MtestArray.append(Mu)
            NmodelArray.append(Nmodel)
            MmodelArray.append(Mmodel)
    NsectionArray = np.array(NsectionArray)
    MsectionArray = np.array(MsectionArray)
    NmodelArray = np.array(NmodelArray)
    MmodelArray = np.array(MmodelArray)
    NtestArray = np.array(NtestArray)
    MtestArray = np.array(MtestArray)
    return NtestArray, NmodelArray, NsectionArray,\
            MtestArray, MmodelArray, MsectionArray


if __name__ == '__main__':
    cfrpdatabase = sio.loadmat('./database/frpdatabase.mat')['cfrpdatabase']
    colArray, Nu, Mu, e0 = testcolarray(cfrpdatabase, failcode=1, fccmodel='tengetal09', code='gb')
    e1testArray, e1modelArray, e1sectionArray,\
        etestArray, emodelArray, esectionArray = test2model_ecc(colArray, Nu, Mu, e0, colmodel='jiangteng13e1', code='gb')
    NtestArray, NmodelArray, NsectionArray,\
            MtestArray, MmodelArray, MsectionArray = test2model_Ru(colArray, Nu, Mu, e0, colmodel='jiangteng13', code='gb')

    sio.savemat('./data/modelerror-cfrp.mat', {'e1testArray': e1testArray,
        'e1modelArray': e1modelArray, 'e1sectionArray': e1sectionArray,
        'NsectionArray': NsectionArray, 'MsectionArray': MsectionArray,
        'NmodelArray': NmodelArray, 'MmodelArray': MmodelArray,
        'NtestArray': NtestArray, 'MtestArray': MtestArray
        })
