import numpy as np

from section import Price, Cost, Material, Geometry
from column import Column

tfsol = np.load('tfsol.npy')
# irrelevant data for column initialization
matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
        'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
price = Price(matprice=matprice)
cost = Cost(price)

def benchmark(tfsol=tfsol, price=price, cost=cost):
    # benchmark data
    h = 400.0
    l = 1000.0
    e0 = 0.15*h
    d2D = 0.80
    fco = 30.0
    ecu = 0.0033
    eco = 0.002
    fr = 400.0
    Er = 200e3
    rhos = 0.01
    ns = 12.
    eh = 1.5/100
    Eh = 230e3
    tft = tfsol[1,0]
    # benchmark column
    # define material
    fc = fco
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, Er=Er, concss=None, reinss=reinss)
    fh = eh*Eh
    mat.setconfinemat(fh, Eh, eh)
    # define geometry
    Ac = np.pi*h**2/4.
    Afb = np.array([0.])    #frp bars
    Astotal = rhos*Ac
    tmp = np.arange(1., ns+1., 2)
    nsportion = np.insert(tmp, [0, ns/2], [0., ns])
    Asportion = (nsportion[1:]-nsportion[:-1])/ns
    As = Astotal*Asportion
    xf = np.array([0.])    #frp bars
    # steel location
    alphas = np.linspace(0., 1., ns/2+1)*np.pi
    rs = d2D*h/2.
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
    column.setcolgeo(l, e0)
    return column

def designcol(tfsol=tfsol, **kwargs):
    benchcol = benchmark()
    # for charateristic design cases
    if kwargs.has_key('e2D') and (kwargs['e2D'] is not None):
        e2D = kwargs['e2D']
        benchcol.e0 = benchcol.geo.h*e2D
    if kwargs.has_key('rhos') and (kwargs['rhos'] is not None):
        rhos = kwargs['rhos']
        rhos0 = np.sum(benchcol.geo.As)/benchcol.geo.Ac
        benchcol.geo.As = benchcol.geo.As*rhos/rhos0
    if kwargs.has_key('d2D') and (kwargs['d2D'] is not None):
        d2D = kwargs['d2D']
        h = benchcol.geo.h
        ns = np.shape(benchcol.geo.xs)[0]*2-2
        alphas = np.linspace(0., 1., ns/2+1)*np.pi
        rs = d2D*h/2.
        xs = h/2. - rs*np.cos(alphas)
        benchcol.geo.xs = xs
    if kwargs.has_key('fco') and (kwargs['fco'] is not None):
        fco = kwargs['fco']
        fcc2fco0 = benchcol.mat.fcc/benchcol.mat.fco
        if np.isclose(fcc2fco0,1.25,atol=1e-3):
            rhoi = 0
        elif np.isclose(fcc2fco0,1.50,atol=1e-3):
            rhoi = 1
        elif np.isclose(fcc2fco0,1.75,atol=1e-3):
            rhoi = 2
        if fco == 30.0:
            fci = 0
        elif fco == 40.0:
            fci = 1
        elif fco == 50.0:
            fci = 2
        tft = tfsol[rhoi,fci]
        benchcol.geo.tft = tft
        benchcol.mat.fc = fco
        benchcol.mat.fco = fco
    if kwargs.has_key('fcc2fco') and (kwargs['fcc2fco'] is not None):
        fcc2fco = kwargs['fcc2fco']
        if fcc2fco == 1.25:
            rhoi = 0
        elif fcc2fco == 1.50:
            rhoi = 1
        elif fcc2fco == 1.75:
            rhoi = 2
        if benchcol.mat.fco == 30.0:
            fci = 0
        elif benchcol.mat.fco == 40.0:
            fci = 1
        elif benchcol.mat.fco == 50.0:
            fci = 2
        tft = tfsol[rhoi,fci]
        benchcol.geo.tft = tft
    benchcol.setfrpcolmat()
    return benchcol


def samplecol(benchcol, **kwargs):
    # for MC simulation
    if kwargs.has_key('h') and (kwargs['h'] is not None):
        hsmp = kwargs['h']
        benchcol.geo.h = hsmp
    if kwargs.has_key('fco') and (kwargs['fco'] is not None):
        fcosmp = kwargs['fco']
        benchcol.mat.fco = fcosmp
        benchcol.mat.fc = fcosmp
    if kwargs.has_key('fs') and (kwargs['fs'] is not None):
        frsmp = kwargs['fs']
        benchcol.mat.fr = frsmp
    if kwargs.has_key('ff') and (kwargs['ff'] is not None):
        ffsmp = kwargs['ff']
        benchcol.mat.eh = ffsmp/benchcol.mat.Eh
    if kwargs.has_key('Ef') and (kwargs['Ef'] is not None):
        eh0 = benchcol.mat.eh
        Eh0 = benchcol.mat.Eh
        Efsmp = kwargs['Ef']
        benchcol.mat.eh = eh0*Eh0/benchcol.mat.Eh
        benchcol.mat.Eh = Efsmp
    benchcol.setfrpcolmat()
    return benchcol

if __name__ == '__main__':
    col0 = benchmark()
    col1 = designcol(fcc2fco=1.75)
    col2 = designcol(fco=40.)
    col3 = designcol(d2D=0.9)
    col4 = designcol(rhos=0.03)
    col5 = designcol(e2D=0.25)
