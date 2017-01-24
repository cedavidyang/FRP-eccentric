import numpy as np
import scipy.integrate as integrate
import sys

RHO_STEEL = 7800.    #kg/m3
RHO_FRP = 1800.
RHO_CONC = 2500.
RHO_SAND = 1600.

class Geometry(object):
    def __init__(self, rtype, h, Ac, Afb, tft, As, xf, xs, bdist):
        """
        rtype: longitudinal reinforcement type, rc for steel rc, frp for frp rc
        Afb: f for frp, b for bar; used for FRP rc
        th: t for thickness, h for hoop; used for FRP-confined concrete
        """
        self.rtype=rtype
        self.h = h
        self.Ac = Ac
        self.Afb = Afb
        self.tft = tft
        self.As = As
        self.xf = xf
        self.xs = xs
        self.bdist = bdist
        if rtype.lower() == 'frp':
            self.Ar = Afb
            self.xr = xf
        elif rtype.lower() == 'rc':
            self.Ar = As
            self.xr = xs
        else:
            print "unknown type, must by frp or rc"
            sys.exit(1)

    def addprop(self,propdict):
        try:
            self.ds = propdict['ds']
        except KeyError:
            pass
        try:
            self.Arepair = propdict['Arepair']
        except KeyError:
            pass
        try:
            self.cover = propdict['cover']
        except KeyError:
            pass

    def addebprop(self, b, h, bf, tf, frpend, shear, span, d, dcmp):
        self.b = b
        self.h = h
        self.bf = bf
        self.tf = tf
        self.frpend = frpend
        self.shear = shear
        self.span = span
        self.d = d
        self.dcmp = dcmp


class Material(object):
    def __init__(self, fc, fr, eco, ecu, eru, Er=None, concss=None, reinss=None):
        self.fc = fc
        self.fr = fr
        self.Er = Er
        self.eco = eco
        self.ecu=ecu
        self.eru=eru
        self.concss = concss
        self.reinss = reinss

    def setconfinemat(self, fh, Eh, eh, hoopss=None):
        self.fh = fh
        self.Eh = Eh
        self.eh = eh
        self.hoopss = hoopss

    def setconfineconcmat(self, fcc, ecu, concss):
        self.fco = self.fc
        self.ecu0 = self.ecu
        self.concss0 = self.concss
        self.fcc = fcc
        self.ecu = ecu
        self.concss = concss


class Cost(object):
    def __init__(self, price=None):
        self.price = price

    def setmatcost(self, conccost, steelcost, fbcost, ftcost):
        self.matcost = conccost+steelcost+fbcost+ftcost

    def settranscost(self, conccost, steelcost, fbcost, ftcost):
        self.transcost = conccost+steelcost+fbcost+ftcost


class Price(object):
    def __init__(self, matprice=None, mntprice=None, failprice=None):
        self.matprice = matprice
        self.mntprice = mntprice
        self.failprice = failprice

    def setmntprice(self, mntprice):
        self.mntprice = mntprice

    def setmatprice(self, matprice):
        self.matprice = matprice

    def setfailprice(self, failprice):
        self.failprice = failprice

    def transprice(self, cargo, distance):
        if cargo.lower() == 'sand':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'water':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'concrete':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'steel':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'fb':
            transprice = 0.035*8*distance
        elif cargo.lower() == 'ft':
            transprice = 0.035*8*distance
        else:
            print "unknow type of cargo"
            sys.exit(1)
        return transprice


def frpdegrade(tday):
    """ t in week"""
    r = (1-0.75)*np.exp(-tday/339.2)+0.75
    return r
