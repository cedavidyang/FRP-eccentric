import numpy as np
import copy

from constants import *
from section import Material, Geometry
from column import Column

def fcc2fco(tf, fco, eco=0.002, ecu=0.0033):
    fr = FRDESIGN
    Er = ERDESIGN
    eh = EPHDESIGN
    Eh = EHDESIGN
    h = HDESIGN
    rhos = RHOSDESIGN
    ns = NSDESIGN
    d2D = DS2DDESIGN
    # define material
    fc = fco
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    # reinss must be a vector function
    reinss = lambda x,Er=Er,fr=fr: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, Er=Er, concss=None, reinss=reinss)
    fh = eh*Eh
    mat.setconfinemat(fh, Eh, eh)
    # define geometry
    Ac = np.pi*h**2/4.
    Afb = np.array([0.])    #frp bars
    tft = tf
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
    fcc = column.mat.fcc
    fco = column.mat.fco
    rho = fcc/fco
    return rho


def benchmark(tfsol, price=price, cost=cost):
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
    reinss = lambda x,Er=Er,fr=fr: np.maximum(np.minimum(x*Er, fr), -fr)
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


def designcol(tfsol, **kwargs):
    benchcol = benchmark(tfsol)
    # for charateristic design cases
    if kwargs.has_key('e2D') and (kwargs['e2D'] is not None):
        e2D = kwargs['e2D']
        benchcol.e0 = benchcol.geo.h*e2D
    if kwargs.has_key('rhos') and (kwargs['rhos'] is not None):
        rhos = kwargs['rhos']
        rhos0 = np.sum(benchcol.geo.As)/benchcol.geo.Ac
        benchcol.geo.As = benchcol.geo.As*rhos/rhos0
        benchcol.geo.Ar = benchcol.geo.As
    if kwargs.has_key('d2D') and (kwargs['d2D'] is not None):
        d2D = kwargs['d2D']
        h = benchcol.geo.h
        ns = np.shape(benchcol.geo.xs)[0]*2-2
        alphas = np.linspace(0., 1., ns/2+1)*np.pi
        rs = d2D*h/2.
        xs = h/2. - rs*np.cos(alphas)
        benchcol.geo.xs = xs
        benchcol.geo.xr = xr
    if kwargs.has_key('fco') and (kwargs['fco'] is not None):
        fco = kwargs['fco']
        rhofcc0 = benchcol.mat.fcc/benchcol.mat.fco
        if np.isclose(rhofcc0,1.25,atol=1e-3):
            rhoi = 0
        elif np.isclose(rhofcc0,1.50,atol=1e-3):
            rhoi = 1
        elif np.isclose(rhofcc0,1.75,atol=1e-3):
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
    if kwargs.has_key('rhofcc') and (kwargs['rhofcc'] is not None):
        rhofcc = kwargs['rhofcc']
        if rhofcc == 1.25:
            rhoi = 0
        elif rhofcc == 1.50:
            rhoi = 1
        elif rhofcc == 1.75:
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
        benchcol.mat.fh = ffsmp
    if kwargs.has_key('Ef') and (kwargs['Ef'] is not None):
        eh0 = benchcol.mat.eh
        Eh0 = benchcol.mat.Eh
        Efsmp = kwargs['Ef']
        benchcol.mat.eh = eh0*Eh0/benchcol.mat.Eh
        benchcol.mat.Eh = Efsmp
    benchcol.setfrpcolmat()
    return benchcol


def loadcol(repcol, code='gb', phif=1.0, gammacc=1.0):
    """To get design Nu,d for loading determination"""
    repcol = copy.deepcopy(repcol)
    fconom = repcol.mat.fc
    frnom = repcol.mat.fr
    fhnom = repcol.mat.fh
    ehnom = repcol.mat.eh
    e0nom = repcol.e0

    if code.lower() == 'gb':
        # When using Teng et al.'s (2009) model, gamma_cc for fcc is equivalent to adding
        # gamma_cc to gamma_c
        repcol.mat.fco = fconom/GAMMA_FC/gammacc
        repcol.mat.fc = fconom/GAMMA_FC/gammacc
        fr = frnom/GAMMA_FS
        Er = repcol.mat.Er
        reinss = lambda x,Er=Er,fr=fr: np.maximum(np.minimum(x*Er, fr), -fr)
        repcol.mat.fr = fr
        repcol.mat.reinss = reinss
        fh = fhnom/GAMMA_FF
        eh = ehnom/GAMMA_FF
        Eh = repcol.mat.Eh
        repcol.mat.setconfinemat(fh, Eh, eh)
        # repcol.e0 = 0.4*e0nom+np.maximum(1.0/30.0*repcol.geo.h, 20.0)

    repcol.setfrpcolmat(model='tengetal09', code=code)

    return repcol


if __name__ == '__main__':
    eco = ECO_GB
    ecu = ECU_GB
    fr = FRDESIGN
    Er = ERDESIGN
    eh = EPHDESIGN
    Eh = EHDESIGN
    h = HDESIGN
    rhos = RHOSDESIGN
    ns = NSDESIGN
    d2D = DS2DDESIGN
    # test of fcc2fco
    fcc2fcoArray = np.array([1.25, 1.50, 1.75])
    fcoArray = np.array([30., 40., 50.])
    tflog = np.empty((3,3),dtype=float)
    for rhoindx, rhoi in enumerate(fcc2fcoArray):
        for fcoindx, fcoi in enumerate(fcoArray):
            tfmin = 0.01*fcoi/eco*h/(2*Eh)
            sol = minimize(lambda x: abs(fcc2fco(x,fcoi,eco,ecu)-rhoi), x0=1.0,
                    bounds=((tfmin,None),),
                    options={'iprint':1})
            tf = sol.x
            tflog[rhoindx, fcoindx] = tf
    # test of design cases
    np.save('tfsol.npy', tflog)
    tfsol = np.load('tfsol.npy')
    col0 = benchmark(tfsol)
    col1 = designcol(tfsol, rhofcc=1.75)
    col2 = designcol(tfsol, fco=40.)
    col3 = designcol(tfsol, d2D=0.9)
    col4 = designcol(tfsol, rhos=0.03)
    col5 = designcol(tfsol, e2D=0.25)
    coldesign = loadcol(col0, code='gb')
