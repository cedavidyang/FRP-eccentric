import numpy as np

from section import Price, Cost, Material, Geometry
from column import Column

if __name__ == '__main__':
    # dummy data
    matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
            'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    price = Price(matprice=matprice)
    cost = Cost(price)

    # unconfined RC column
    # define material
    fc = 20.1
    fr = 335.
    Er = 200e3
    ecu = 0.0033
    eco = 0.002
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    def concss(ec, eco=eco, fc=fc):    # must be scalar function
        if ec<=eco:
            sc = fc*(1.-(1.-ec/eco)**(2-20./60))
        else:
            sc = fc
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 600.
    Ac = np.pi*h**2/4.
    Afb = np.array([0.])
    tft = 0.
    Astotal = 0.02*Ac
    As = np.array([Astotal/12., Astotal/6., Astotal/6., Astotal/6.,
        Astotal/6., Astotal/6., Astotal/12.])
    xf = np.array([0.])
    xs = np.array([50, 83.5, 175, 300, 425, 516.5, 550])
    def bdist(x, h=h, b=h):
        if x<0 or x>h:
            bsec = 0.
        else:
            bsec = 2*np.sqrt((h/2)**2-(x-h/2)**2)
        return bsec
    geo = Geometry('rc', h, Ac, Afb, tft, As, xf, xs, bdist)
    # define beam
    column0 = Column(geo=geo, mat=mat, cost=cost)
    N0 = column0.axial_compression()
    print 'axial compression = {} kN'.format(N0/1e3)
    dummy,M0,cbending = column0.capacity(0)
    print 'pure bending = {} kNm'.format(M0/1e6)


    # verification 1 (Jiang 2008 Fig. 5.12(b))
    # define material
    fc = 20.1
    fr = 335.
    Er = 200e3
    ecu = 0.0033
    eco = 0.002
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss=None, reinss=reinss)
    Eh = 230e3; eh = 0.0075; fh = Eh*eh
    mat.setconfinemat(fh, Eh, eh)
    # define geometry
    h = 600.
    Ac = np.pi*h**2/4.
    Afb = np.array([0.])
    tft = 0.63
    Astotal = 0.02*Ac
    As = np.array([Astotal/12., Astotal/6., Astotal/6., Astotal/6.,
        Astotal/6., Astotal/6., Astotal/12.])
    xf = np.array([0.])
    xs = np.array([50, 83.5, 175, 300, 425, 516.5, 550])
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
    print('fcc/fco = {}'.format(column.mat.fcc/column.mat.fco))

    # section analysis
    Ncc0 = column.axial_compression(confine='yes')
    NarraySection = np.arange(1e-3, 1-1e-3, 0.1)*Ncc0
    MarraySection = []
    ciarray = []
    for Ni in NarraySection:
        Nsec,Msec,ci = column.capacity(Ni)
        MarraySection.append(Msec)
        ciarray.append(ci)
    MarraySection = np.array(MarraySection)

    # jiang and teng (2013) model
    thetaArray = np.linspace(0.1, 0.75, 50)
    NarrayModel = []
    MarrayModel = []
    for thetai in thetaArray:
        Nmodeli = column.modelNmat(thetai)
        Mmodeli = column.modelMmat(thetai)
        NarrayModel.append(Nmodeli)
        MarrayModel.append(Mmodeli)
    NarrayModel = np.array(NarrayModel)
    MarrayModel = np.array(MarrayModel)
    indx =  (NarrayModel<0) | (MarrayModel<0)
    NarrayModel = NarrayModel[~indx]
    MarrayModel = MarrayModel[~indx]

    import scipy.io as sio
    sio.savemat('./data/jiangteng2013NM.mat', {'N0': N0, 'M0':M0,
        'NarraySection':NarraySection, 'MarraySection':MarraySection,
        'NarrayModel':NarrayModel, 'MarrayModel':MarrayModel})
