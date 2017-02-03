import numpy as np
import scipy.io as sio
from scipy.optimize import minimize

from section import Price, Cost, Material, Geometry
from column import Column
from constants import *

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
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
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

if __name__ == '__main__':
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
    np.save('tfsol.npy', tflog)
