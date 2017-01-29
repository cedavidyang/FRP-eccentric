import numpy as np
import scipy.io as sio
from scipy.optimize import minimize
from scipy.interpolate import interp1d

from section import Price, Cost, Material, Geometry
from column import Column

hcol = 1
l2hcol = 2
tftcol = 3
Ehcol = 4
ehcoupon = 6
frcol = 11
Ercol = 12
rhoscol = 13
covercol = 15
nlcol = 17
fccol = 18
e0col = 19
Nucol = 20
Mucol = 21
ehrup = 23
fmcol = 25


if __name__ == '__main__':
    cfrpdatabase = sio.loadmat('./database/frpdatabase.mat')['cfrpdatabase']
    e1modelArray = []
    e1sectionArray = []
    e1section0Array = []
    e1testArray = []
    failmode = np.array([columnData[fmcol] for columnData in cfrpdatabase])
    for icol, columnData in enumerate(cfrpdatabase[failmode==1]):
        # irrelevant data for column initialization
        matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
                'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
        price = Price(matprice=matprice)
        cost = Cost(price)
        # define material
        fc = columnData[fccol]
        fr = columnData[frcol]
        Er = columnData[Ercol]*1e3
        ecu = 0.0033
        eco = 0.002
        # rupture strain of steel is assmued to be 100 times larger than yield strain
        eru = fr/Er*100
        # reinss must be a vector function
        reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
        mat = Material(fc, fr, eco, ecu, eru, Er=Er, concss=None, reinss=reinss)
        eh = columnData[ehcoupon]
        Eh = columnData[Ehcol]*1e3
        fh = eh*Eh
        mat.setconfinemat(fh, Eh, eh)
        # define geometry
        h = columnData[hcol]
        Ac = np.pi*h**2/4.
        Afb = np.array([0.])    #frp bars
        tft = columnData[tftcol]
        Astotal = columnData[rhoscol]*Ac
        ns = columnData[nlcol]
        tmp = np.arange(1., ns+1., 2)
        nsportion = np.insert(tmp, [0, ns/2], [0., ns])
        Asportion = (nsportion[1:]-nsportion[:-1])/ns
        As = Astotal*Asportion
        xf = np.array([0.])    #frp bars
        # steel location
        alphas = np.linspace(0., 1., ns/2+1)*np.pi
        rs = h/2. - columnData[covercol]
        xs = h/2. - rs*np.cos(alphas)
        def bdist(x, h=h, b=h):
            if x<0 or x>h:
                bsec = 0.
            else:
                bsec = 2*np.sqrt((h/2)**2-(x-h/2)**2)
            return bsec
        geo = Geometry('rc', h, Ac, Afb, tft, As, xf, xs, bdist)
        # define beam
        column = Column(geo=geo, mat=mat, cost=cost)
        column.setfrpcolmat()
        Nu = columnData[Nucol]*1e3
        Mu = columnData[Mucol]*1e6
        e0 = columnData[e0col]
        etotal = Mu/Nu
        e1test = etotal - e0
        # extra eccentricity based on Jiang and Teng (2013)
        l = columnData[l2hcol]*columnData[hcol]
        column.setcolgeo(l, e0)
        Ncolmodel, Mcolmodel, e1model = column.colcapacitymodel()
        print 'model: {} out of {}'.format(icol+1, np.shape(cfrpdatabase[failmode==1])[0])
        Ncolsection, Mcolsection, e1section = column.colcapacitysection_fittedcurve()
        # add to log
        e1testArray.append(e1test)
        e1modelArray.append(e1model)
        e1sectionArray.append(e1section)
        print 'section: {} out of {}'.format(icol+1, np.shape(cfrpdatabase[failmode==1])[0])

    e1testArray = np.array(e1testArray)
    e1modelArray = np.array(e1modelArray)
    e1sectionArray = np.array(e1sectionArray)


    MsectionArray = []
    NsectionArray = []
    MmodelArray = []
    NmodelArray = []
    MtestArray = []
    NtestArray = []
    for icol,columnData in enumerate(cfrpdatabase[failmode==1]):
        # irrelevant data for column initialization
        matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
                'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
        price = Price(matprice=matprice)
        cost = Cost(price)
        # define material
        fc = columnData[fccol]
        fr = columnData[frcol]
        Er = columnData[Ercol]*1e3
        ecu = 0.0033
        eco = 0.002
        # rupture strain of steel is assmued to be 100 times larger than yield strain
        eru = fr/Er*100
        # reinss must be a vector function
        reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
        mat = Material(fc, fr, eco, ecu, eru, Er=Er, concss=None, reinss=reinss)
        eh = columnData[ehcoupon]
        Eh = columnData[Ehcol]*1e3
        fh = eh*Eh
        mat.setconfinemat(fh, Eh, eh)
        # define geometry
        h = columnData[hcol]
        Ac = np.pi*h**2/4.
        Afb = np.array([0.])    #frp bars
        tft = columnData[tftcol]
        Astotal = columnData[rhoscol]*Ac
        ns = columnData[nlcol]
        tmp = np.arange(1., ns+1., 2)
        nsportion = np.insert(tmp, [0, ns/2], [0., ns])
        Asportion = (nsportion[1:]-nsportion[:-1])/ns
        As = Astotal*Asportion
        xf = np.array([0.])    #frp bars
        # steel location
        alphas = np.linspace(0., 1., ns/2+1)*np.pi
        rs = h/2. - columnData[covercol]
        xs = h/2. - rs*np.cos(alphas)
        def bdist(x, h=h, b=h):
            if x<0 or x>h:
                bsec = 0.
            else:
                bsec = 2*np.sqrt((h/2)**2-(x-h/2)**2)
            return bsec
        geo = Geometry('rc', h, Ac, Afb, tft, As, xf, xs, bdist)
        # define beam
        column = Column(geo=geo, mat=mat, cost=cost)
        column.setfrpcolmat()
        Nu = columnData[Nucol]*1e3
        Mu = columnData[Mucol]*1e6
        e0 = columnData[e0col]
        etotal = Mu/Nu
        e1test = etotal - e0
        # section analysis
        Nccu = column.axial_compression(confine='yes')
        NarraySection = np.linspace(1e-3, 1-1e-3, 50)*Nccu
        MarraySection = []
        ciarray = []
        for Ni in NarraySection:
            Nsec,Msec,ci = column.capacity(Ni)
            MarraySection.append(Msec)
            ciarray.append(ci[0])
        MarraySection = np.array(MarraySection)
        fNM = interp1d(NarraySection, MarraySection, bounds_error=False,
                fill_value=(MarraySection[0], MarraySection[-1]))
        # def desection(Niscale, scale, column=column, etotal=etotal):
            # Ni = Niscale*scale
            # Nsec,Msec,ci = column.capacity(Ni)
            # return abs(Msec/Nsec-etotal)
        def desection(Niscale, scale, fNM=fNM, etotal=etotal):
            Ni = Niscale*scale
            Msec = fNM(Ni)
            return abs(Msec/Ni-etotal)
        scale = 1e3
        sol1 = minimize(desection, Nu/scale, args=(scale,),
                bounds=((0,Nccu/scale),), options={'iprint':-1})
        Nsection = sol1.x[0]*scale
        Msection = Nsection*etotal
        print 'section: {} out of {}'.format(icol+1, np.shape(cfrpdatabase[failmode==1])[0])
        # jiang and teng (2013) model
        def demodel(thetai, column=column, etotal=etotal):
            Nmodel = column.modelNmat(thetai)
            Mmodel = column.modelMmat(thetai)
            return abs(Mmodel/Nmodel-etotal)
        sol2 = minimize(demodel, (0.1+0.75)/2., bounds=((0.1,0.75),))
        theta = sol2.x[0]
        Nmodel = column.modelNmat(theta)
        Mmodel = column.modelMmat(theta)
        print 'model: {} out of {}'.format(icol+1, np.shape(cfrpdatabase[failmode==1])[0])
        # add to log
        NtestArray.append(Nu)
        MtestArray.append(Mu)
        NsectionArray.append(Nsection)
        MsectionArray.append(Msection)
        NmodelArray.append(Nmodel)
        MmodelArray.append(Mmodel)
    NsectionArray = np.array(NsectionArray)
    MsectionArray = np.array(MsectionArray)
    NmodelArray = np.array(NmodelArray)
    MmodelArray = np.array(MmodelArray)
    NtestArray = np.array(NtestArray)
    MtestArray = np.array(MtestArray)

    sio.savemat('./data/modelerror.mat', {'e1testArray': e1testArray,
        'e1modelArray': e1modelArray, 'e1sectionArray': e1sectionArray,
        'NsectionArray': NsectionArray, 'MsectionArray': MsectionArray,
        'NmodelArray': NmodelArray, 'MmodelArray': MmodelArray,
        'NtestArray': NtestArray, 'MtestArray': MtestArray
        })
