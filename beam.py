import numpy as np
import scipy.integrate as integrate

from section import Price, Cost, Material, Geometry
from section import RHO_FRP, RHO_STEEL, RHO_CONC
from section import frpdegrade

import sys
import copy

class Beam(object):
    def __init__(self, geo=None, mat=None, cost=None):
        self.geo=geo
        self.mat=mat
        self.cost = cost

    def setgeo(self, geo):
        self.geo=geo

    def setmat(self, mat):
        self.mat=mat

    def setcost(self, cost):
        self.cost=cost

    def matcost(self):
        if self.geo is None or self.mat is None:
            print "geometry and material are not defined"
            sys.exit(1)
        if self.cost.price is None:
            print "Cost is not initiated properly"
            sys.exit(1)
        pconc = self.cost.price.matprice['Cconc']
        psteel = self.cost.price.matprice['Csteel']
        pfb = self.cost.price.matprice['Cfb']
        pft = self.cost.price.matprice['Cft']
        conccost = self.geo.Ac*1e-6*pconc
        steelcost = np.sum(self.geo.As)*1e-6*RHO_STEEL*1e-3*psteel
        fbcost = np.sum(self.geo.Afb)*1e-6*RHO_FRP*1e-3*pfb
        ftcost = self.geo.Aft*1e-6*RHO_FRP*1e-3*pft
        self.cost.setmatcost(conccost, steelcost, fbcost, ftcost)
        return self.cost.matcost

    def transcost(self, distance):
        if self.geo is None or self.mat is None:
            print "geometry and material are not defined"
            sys.exit(1)
        if self.cost.price is None:
            print "Cost is not initiated properly"
            sys.exit(1)
        if self.geo.rtype.lower() == 'rc':
            concmass = self.geo.Ac*1e-6*RHO_CONC*1e-3
            # sandratio = 250./903
            # waterratio = 1./301
            # sandcost = sandratio*concmass*self.cost.price.transprice('sand', distance)
            # watercost = waterratio*concmass*self.cost.price.transprice('water', distance)
            # conccost = sandcost+watercost
            conccost = concmass*self.cost.price.transprice('concrete', distance)
            steelmass = np.sum(self.geo.As)*1e-6*RHO_STEEL*1e-3
            steelcost = steelmass*self.cost.price.transprice('steel', distance)
            fbcost = 0.
            ftcost = 0.
        else:
            sandratio = 250./903
            waterratio = 1./301
            concmass = self.geo.Ac*1e-6*RHO_CONC*1e-3
            conccost = concmass*(1-sandratio-waterratio)*self.cost.price.transprice('concrete', distance)
            steelcost = 0.
            fbmass = np.sum(self.geo.Afb)*1e-6*RHO_FRP*1e-3
            fbcost = fbmass*self.cost.price.transprice('fb', distance)
            ftmass = self.geo.Aft*1e-6*RHO_FRP*1e-3
            ftcost = ftmass*self.cost.price.transprice('ft', distance)
        self.cost.settranscost(conccost, steelcost, fbcost, ftcost)
        return self.cost.transcost

    def directmntcost(self, mntplan):
        pdirect = self.cost.price.mntprice['direct']
        pindirect = self.cost.price.mntprice['indirect']
        v = self.cost.price.mntprice['discount']
        nmnt = mntplan['Mt'].size

        # price data
        Cwb = pdirect['Cwb']
        Crm = pdirect['Crm']
        Cbc = pdirect['Cwb']
        Ctp = pdirect['Crm']
        # geo data
        Arm2 = self.geo.Arepair*1e-6    # in m2
        cover = self.geo.cover    # in mm
        Vcover = Arm2*(cover/10)
        # directcost
        directcost = 0.
        for t in mntplan['Mt']:
            cost = ((Cwb+Crm)*Vcover + (Cbc+Ctp)*Arm2)/(1+v)**t
            directcost += cost
        self.cost.mntcost = directcost
        return self.cost.mntcost

    def getmoment(self, nd=100, carray=None):
        """ nh: number of discretiztions along beam (effective) depth d
            carray: given carray
        """
        xr = self.geo.xr
        Ar = self.geo.Ar
        dmin = np.min(xr)
        dmax = np.max(xr)
        bdist = self.geo.bdist    # distribution of concrete width
        ecu = self.mat.ecu
        eru = self.mat.eru
        concss=self.mat.concss    #concrete stress-strain relation
        reinss = self.mat.reinss    # reinforcement stress-strain relation

        if carray is None:
            carray = np.linspace(dmin/nd, dmin, num=nd)
        rfmin = 1.e16
        for c in carray:
            # assume concrete crush, determine steel/frpbar strain
            er = (dmax-c)/c*ecu
            if er>eru:    # control by reinforcement failure
                ec = c/(dmax-c)*eru
            else:
                ec = ecu
            ecdist = lambda x: -ec/c*x+ec
            erdist = lambda x: (x-c)/c*ec
            Fc,err = integrate.quad(lambda x: concss(ecdist(x))*bdist(x), 0, c)
            Fr = np.sum(reinss(erdist(xr))*Ar)
            rf = np.abs(Fc-Fr)
            if rf<rfmin:
                rfmin = rf
                csol = c
                ecsol = ec
                ecdistsol = lambda x: -ecsol/csol*x+ecsol
                erdistsol = lambda x: (x-csol)/csol*ecsol

        # calculate the capacity
        Mc,err = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x)*x, 0, csol)
        Mr = np.sum(reinss(erdistsol(xr))*Ar*xr)
        M = np.abs(Mc-Mr)
        self.capacity = M
        return M,csol

    def lifecycleR(self, life, mntplan):
        ds = self.geo.ds    # rebar diameter
        tarray = np.arange(1, life+1, dtype=float)
        if self.geo.rtype.lower() == 'rc':
            darray = ds*np.ones(tarray.shape)
            Rarray = np.ones(tarray.shape)
            tinit = 10.
            rcorr = 0.127    # mm/yr
            for i,t in enumerate(tarray):
                if t<tinit:
                    darray[i] = ds
                    Rarray[i] = 1.0
                else:
                    # look for the previous maintenance
                    if mntplan['Mt'] is None:
                        indx=0
                    else:
                        indx = np.searchsorted(mntplan['Mt'], t)
                    if indx==0:    # no maintenance is applied
                        dst = ds-rcorr*(t-tinit)
                        darray[i] = dst
                        Rarray[i] = (dst/ds)**2
                    else:
                        tlastmnt = mntplan['Mt'][indx-1]
                        dt = t-tlastmnt
                        if dt<mntplan['Me']:
                            darray[i] = darray[i-1]
                            Rarray[i] = Rarray[i-1]
                        else:
                            dst = darray[i-1]-rcorr
                            darray[i] = dst
                            Rarray[i] = (dst/ds)**2
        elif self.geo.rtype.lower() == 'frp':
            tday = tarray*360
            Rarray = frpdegrade(tday)
        return Rarray


    def lifecycleR_exact(self, life, mntplan):
        ds = self.geo.ds    # rebar diameter
        tarray = np.arange(1, life+1, dtype=float)
        beamclone = copy.deepcopy(self)
        R0 = self.getmoment()[0]
        fr0 = self.mat.fr
        As0 = self.geo.As
        if self.geo.rtype.lower() == 'rc':
            darray = ds*np.ones(tarray.shape)
            Rarray = np.ones(tarray.shape)
            tinit = 10.
            rcorr = 0.127    # mm/yr
            for i,t in enumerate(tarray):
                if t<tinit:
                    darray[i] = ds
                    Rarray[i] = 1.0
                else:
                    # look for the previous maintenance
                    if mntplan['Mt'] is None:
                        indx=0
                    else:
                        indx = np.searchsorted(mntplan['Mt'], t)
                    if indx==0:    # no maintenance is applied
                        dst = ds-rcorr*(t-tinit)
                        darray[i] = dst
                        # Rarray[i] = (dst/ds)**2
                        beamclone.geo.Ar = As0*(dst/ds)**2
                        Rt = beamclone.getmoment()[0]
                        Rarray[i] = Rt/R0
                    else:
                        tlastmnt = mntplan['Mt'][indx-1]
                        dt = t-tlastmnt
                        if dt<mntplan['Me']:
                            darray[i] = darray[i-1]
                            Rarray[i] = Rarray[i-1]
                        else:
                            dst = darray[i-1]-rcorr
                            darray[i] = dst
                            # Rarray[i] = (dst/ds)**2
                            beamclone.geo.Ar = As0*(dst/ds)**2
                            Rt = beamclone.getmoment()[0]
                            Rarray[i] = Rt/R0
        elif self.geo.rtype.lower() == 'frp':
            tday = tarray*360
            # Rarray = frpdegrade(tday)
            Rarray = np.ones(tarray.shape)
            rofarray = frpdegrade(tday)
            for i,ro in enumerate(rofarray):
                beamclone.mat.fr = fr0*ro
                Rt = beamclone.getmoment()[0]
                Rarray[i] = Rt/R0
        return Rarray

def getrcbeam(price):
    cost = Cost(price)
    fc = 30.
    ecu = 0.003
    fr = 420.; Er = 200000.
    eru = fr/Er*100.
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.003-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=800.
    Ac = 400.*(800.-200)+2500.*200.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([4*819.+2*509])
    xf = np.array([0.])
    xs = np.array([700])
    def bdist(x):
        if x<0 or x>800:
            b = 0.
        elif x>=0 and x<200:
            b = 2500.
        else:
            b = 400.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    geo.addprop({'ds':25.4})
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam


def getrcbeam_cice2016(price=None):
    cost = Cost(price)
    fc = 30.
    ecu = 0.003
    fr = 420.; Er = 200000.
    eru = fr/Er*100.
    Ec = 4700*np.sqrt(fc)
    fcp = 0.85*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.003-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=790.
    Ac = 400.*(790.-190)+2600.*190.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([4*819.+2*819.])
    xf = np.array([0.])
    xs = np.array([700])
    def bdist(x):
        if x<0 or x>790.:
            b = 0.
        elif x>=0 and x<190.:
            b = 2600.
        else:
            b = 400.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    geo.addprop({'ds':32.26})
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam


def getfrpbeam(price):
    cost = Cost(price)
    fc = 30.
    ecu = 0.003
    fr = 550.*0.2; Er = 45000.
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.003-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=800.
    Ac = 400.*(800.-200)+2500.*200.
    Afb = np.array([4*819.+4*661.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([720.])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>800:
            b = 0.
        elif x>=0 and x<200:
            b = 2500.
        else:
            b = 400.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    geo.addprop({'ds':32.26})
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam


def getfrpbeam_cice2016(price=None):
    cost = Cost(price)
    fc = 30.
    ecu = 0.003
    fr = 550.; Er = 45000.
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.003-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=790.
    Ac = 400.*(790.-190.)+2600.*190.
    # Afb = np.array([4*819.+4*661.])
    Afb = np.array([8*819.+0*645.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([710.])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>790:
            b = 0.
        elif x>=0 and x<190.:
            b = 2600.
        else:
            b = 400.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    geo.addprop({'ds':32.26})
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam


def getfrpcreep_cice2016(price=None):
    cost = Cost(price)
    fc = 30.
    ecu = 0.003
    fr = 550.*0.2; Er = 45000.
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.003-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=790.
    Ac = 400.*(790.-190)+2600.*190.
    Afb = np.array([8*819.+0*645.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([710.])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>790:
            b = 0.
        elif x>=0 and x<190:
            b = 2600.
        else:
            b = 400.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    geo.addprop({'ds':32.26})
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    return beam

if __name__ == '__main__':
    matprice = {'Cconc': 104.57*9, 'Csteel': 1.16e3*9, 'Cfb': 1.16e3*9, 'Cft': 1.16e3*9}
    price = Price(matprice=matprice)
    cost = Cost(price)
    ## RC example ex4.1, 4.2, 4.10, 4.11, 4.12 (Shu Shilin concrete textbook) checked
    # define material
    fc = 11.9
    fr = 300; Er = 200e3
    ecu = 0.0033
    eco = 0.002
    eru = fr/Er*100
    def concss(ec):
        if ec<=eco:
            sc = fc*(1.-(1.-ec/eco)**(2-25./50))
        else:
            sc = fc
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=500.
    Ac = 200.*(500.-120)+500.*120.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([2454.])
    xf = np.array([0.])
    xs = np.array([442.5])
    def bdist(x):
        if x<0 or x>500:
            b = 0.
        elif x>=0 and x<120:
            b = 500.
        else:
            b = 200.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    print 'material cost = {}'.format(beam.matcost())
    M,c = beam.getmoment()
    print 'capacity = {} kN-m; neutral axis at x={}'.format(M/1e6, c)


    ## FRP example example beam in ACI440.1R-06 checked (concrete stress-strain curve
    ## in Wight and MacGregor Textbook has been used)
    # define material
    fc = 27.6
    fr = 496.; Er = 44800.
    ecu = 0.003
    eru = fr/Er
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    reinss = lambda x: np.minimum(x*Er, fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h=305.
    Ac = 305.*178.
    Afb = np.array([400.])
    Aft = 0.
    As = np.array([0.])
    xf = np.array([250.])
    xs = np.array([0.])
    def bdist(x):
        if x<0 or x>305:
            b = 0.
        else:
            b = 178.
        return b
    geo = Geometry('frp', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    beam = Beam(geo=geo, mat=mat, cost=cost)
    print 'material cost = {}'.format(beam.matcost())
    M,c = beam.getmoment()
    print 'capacity = {} kN-m; neutral axis at x={}'.format(M/1e6, c)
