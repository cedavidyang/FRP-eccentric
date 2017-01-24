import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
from scipy.optimize import root
from scipy.interpolate import interp1d

from section import Price, Cost, Material, Geometry
from beam import Beam
import sys

class Column(Beam):
    def setcolgeo(self, l, e0, colend='pinned'):
        self.l = l
        self.e0 = e0
        self.colend = colend

    def setfrpcolmat(self, model='tengetal09'):
        if model.lower() == 'tengetal09':
            Esec0 = self.mat.fc/self.mat.eco
            rhoK = 2*self.mat.Eh*self.geo.tft/(Esec0*self.geo.h)
            rhoE = self.mat.eh/self.mat.eco
            if rhoK<0.01:
                fcc2fco = 1.
            else:
                fcc2fco = 1+3.5*(rhoK-0.01)*rhoE
            fco = self.mat.fc
            fcc = fcc2fco * fco
            # 1.75 is replaced by 1.65 for GB
            ecu2eco = 1.65+6.5*rhoK**0.8*rhoE**1.45
            ecu = ecu2eco*self.mat.eco
            Ec = 1000*fco
            E2 = (fcc-fco)/ecu
            et = 2*fco/(Ec-E2)
            def concss(ec, fco=fco, fcc=fcc, E2=E2, Ec=Ec, et=et, ecu=ecu):    # must be scalar function
                if ec>=0 and ec<=et:
                    sc = Ec*ec - (Ec-E2)**2/(4*fco)*ec**2
                elif ec>et and ec<ecu:
                    sc = fco + E2*ec
                else:
                    sc = 0.
                return sc
            self.mat.setconfineconcmat(fcc, ecu, concss)
        elif model.lower() == 'lamteng03':
            print "available in future version"
            sys.exit(1)

    def axial_compression(self, confine='no'):
        if confine is 'yes':
            if self.geo.rtype.lower() == 'rc':
                # consider reinforcement contribution
                Nr = np.sum(self.geo.Ar*self.mat.fr)
                Nc = (self.geo.Ac-np.sum(self.geo.Ar))*self.mat.fcc
            else:
                Nr = 0.0
                Nc = self.geo.Ac*self.mat.fcc
        else:
            if self.geo.rtype.lower() == 'rc':
                # consider reinforcement contribution
                Nr = np.sum(self.geo.Ar*self.mat.fr)
                Nc = (self.geo.Ac-np.sum(self.geo.Ar))*self.mat.fc
            else:
                Nr = 0.0
                Nc = self.geo.Ac*self.mat.fc
        return Nr+Nc

    def capacity(self, N):
        """ nh: number of discretiztions along beam (effective) depth d"""
        h = self.geo.h
        xr = self.geo.xr
        Ar = self.geo.Ar
        d = np.max(xr)
        bdist = self.geo.bdist    # distribution of concrete width
        ecu = self.mat.ecu
        eru = self.mat.eru
        if self.geo.rtype.lower() == 'rc':
            negeru = -eru
        else:
            negeru = -eru
        concss=self.mat.concss    #concrete stress-strain relation
        concss_vector = np.vectorize(concss)
        reinss = self.mat.reinss    # reinforcement stress-strain relation

        rfmin = 1.e16
        def objfunc(c):
            # assume concrete crush, determine steel/frpbar strain
            er = (d-c)/c*ecu
            if er>eru or er<negeru:    # control by reinforcement failure
                ec = c/(d-c)*eru
            else:
                ec = ecu
            ecdist = lambda x: -ec/c*x+ec
            erdist = lambda x: (x-c)/c*ec
            Fc,err = integrate.quad(lambda x: concss(ecdist(x))*bdist(x), 0, c)
            Fccancel = np.sum(concss_vector(ecdist(xr))*Ar)
            Fr = np.sum(reinss(erdist(xr))*Ar)
            rf = np.abs( (Fc-Fccancel)-Fr-N)
            return rf
        c0 = 0.5*d
        res = minimize(objfunc, c0, method='L-BFGS-B', bounds=((0.1*d,None),))
        # res = minimize(objfunc, c0, bounds=((0.1*d,None),))
        csol = res.x

        er = (d-csol)/csol*ecu
        if er>eru or er<negeru:    # control by reinforcement failure
            ecsol = csol/(d-csol)*eru
        else:
            ecsol = ecu
        ecdistsol = lambda x: -ecsol/csol*x+ecsol
        erdistsol = lambda x: (x-csol)/csol*ecsol

        # calculate the capacity
        Fc,errN = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x), 0, csol)
        Fccancel = np.sum(concss_vector(ecdistsol(xr))*Ar)
        Fr = np.sum(reinss(erdistsol(xr))*Ar)
        Ncomp = np.abs( (Fc-Fccancel)-Fr )
        Mc,errM = integrate.quad(lambda x: concss(ecdistsol(x))*bdist(x)*x, 0, csol)
        Mccancel = np.sum(concss_vector(ecdistsol(xr))*Ar*xr)
        Mr = np.sum(reinss(erdistsol(xr))*Ar*xr)
        M = np.abs( (Mc-Mccancel)-Mr-Ncomp*h/2. )
        return Ncomp,M,csol

    def modelNmat(self, theta, model='jiangteng13'):
        """ theta: polar coordinate of neutral axis in circular columns
        0.1 \leq thea \leq 0.75 (Jiang & Teng 2013 Eq. 22d and 22e)
        N expression is based on Jiang and Teng (2013), Eq. 22a
        """
        fcc = self.mat.fcc
        fco = self.mat.fco
        A = np.pi*self.geo.h**2/4.
        fy = self.mat.fr
        As = np.sum(self.geo.As)
        alpha1 = 1.17-0.2*fcc/fco
        thetac = 1.25*theta-0.125
        if thetac<0:
            thetac = 0.
        elif thetac>1:
            thetac = 1.
        thetat = 1.125-1.5*theta
        if thetat<0:
            thetat = 0.
        elif thetat>1:
            thetat = 1.
        # approximation of N from section analysis
        Nmat = theta*alpha1*fcc*A*(1-np.sin(2*np.pi*theta)/(2*np.pi*theta)) + \
                (thetac-thetat)*fy*As
        return Nmat

    def modelMmat(self, theta, model='jiangteng13'):
        """ theta: polar coordinate of neutral axis in circular columns
        0.1 \leq thea \leq 0.75 (Jiang & Teng 2013 Eq. 22d and 22e)
        N expression is based on Jiang and Teng (2013), Eq. 22b
        """
        fcc = self.mat.fcc
        fco = self.mat.fco
        A = np.pi*self.geo.h**2/4.
        R = self.geo.h/2.
        Rs = self.geo.h/2.-self.geo.xr[0]
        fy = self.mat.fr
        As = np.sum(self.geo.As)
        alpha1 = 1.17-0.2*fcc/fco
        thetac = 1.25*theta-0.125
        if thetac<0:
            thetac = 0.
        elif thetac>1:
            thetac = 1.
        thetat = 1.125 - 1.5*theta
        if thetat<0:
            thetat = 0.
        elif thetat>1:
            thetat = 1.
        Mmat = 2./3.*alpha1*fcc*A*R*np.sin(np.pi*theta)**3/np.pi + \
                fy*As*Rs*(np.sin(np.pi*thetac)+np.sin(np.pi*thetat))/np.pi
        return Mmat

    def colcapacitymodel(self, model='jiangteng13'):
        ecu = self.mat.ecu
        ey = self.mat.fr/self.mat.Er
        D = self.geo.h
        d = self.geo.h - 2*self.geo.xr[0]
        phibal = 2*(ecu+ey)/(D+d)
        fcc = self.mat.fcc
        Ac = self.geo.Ac
        Nbal = 0.8*fcc*Ac
        Nccu = self.axial_compression(confine='yes')
        solthetalb = minimize(lambda x: abs(self.modelNmat(x)-0.1*Nccu),
                (0.1+0.75)/2., bounds=((0.1,0.75),))
        thetalb = solthetalb.x[0]
        soltheta0 = minimize(lambda x: abs(self.modelNmat(x)-Nbal),
                (0.1+0.75)/2., bounds=((0.1,0.75),))
        theta0 = soltheta0.x[0]
        def solvecolumn(theta, column=self, Nbal=Nbal, phibal=phibal):
            e0 = column.e0
            l = column.l
            Nu = column.modelNmat(theta)
            xi1 = Nbal/Nu
            if xi1>1: xi1=1
            Mu1 = Nu*(e0+l**2/np.pi**2*xi1*phibal)
            Mu2 = column.modelMmat(theta)
            return abs(Mu1-Mu2)
        x0candidate = np.linspace(thetalb, 0.75, 50)
        dM = np.array([solvecolumn(thetai) for thetai in x0candidate])
        x0 = x0candidate[np.argmin(dM)]
        sol = minimize(solvecolumn, x0, bounds=((thetalb, 0.75),))
        theta = sol.x[0]
        Ncol = self.modelNmat(theta)
        Mcol = self.modelMmat(theta)
        e1 = Mcol/Ncol - self.e0
        return Ncol, Mcol, e1

    def colcapacitysection(self):
        Nccu = self.axial_compression(confine='yes')
        def solvecolumn(N, column=self):
            e0 = column.e0
            l = column.l
            d = np.max(column.geo.xr)
            ecu = column.mat.ecu
            eru = column.mat.eru
            Nsec, Msec, csol = column.capacity(N)
            csol = csol[0]
            er = (d-csol)/csol*ecu
            eru = column.mat.eru
            if er>eru or er<-eru:    # control by reinforcement failure
                ecf = csol/(d-csol)*eru
            else:
                ecf = ecu
            phiu = ecf/csol
            Mu1 = N*(e0+l**2/np.pi**2*phiu)
            Mu2 = Msec
            return abs(Mu1-Mu2)
        x0candidate = np.linspace(1e-2, 1-1e-2, 10)*Nccu
        dM = np.array([solvecolumn(Nui) for Nui in x0candidate])
        x0 = x0candidate[np.argmin(dM)]
        sol = minimize(solvecolumn, x0, bounds=((0.01*Nccu, Nccu),))
        Ncol = sol.x[0]
        Ncol, Mcol, csol = self.capacity(Ncol)
        e1 = Mcol/Ncol - self.e0
        return Ncol, Mcol, e1

    def colcapacitysection_fittedcurves(self, NM=None, Nc=None):
        Nccu = self.axial_compression(confine='yes')
        if NM is None or Nc is None:
            NarraySection = np.linspace(1e-3, 1-1e-3, 50)*Nccu
            MarraySection = []
            ciarray = []
            for Ni in NarraySection:
                Nsec,Msec,ci = self.capacity(Ni)
                MarraySection.append(Msec)
                ciarray.append(ci[0])
            MarraySection = np.array(MarraySection)
            fNM = interp1d(NarraySection, MarraySection, bounds_error=False,
                    fill_value=(MarraySection[0], MarraySection[-1]))
            fNci = interp1d(NarraySection, ciarray, bounds_error=False,
                    fill_value=(ciarray[0], ciarray[-1]))
        else:
            fNM = NM
            fNci = Nc
        def solvecolumn(N, column=self, fNM=fNM, fNci=fNci):
            e0 = column.e0
            l = column.l
            d = np.max(column.geo.xr)
            ecu = column.mat.ecu
            eru = column.mat.eru
            Msec = fNM(N)
            csol = fNci(N)
            er = (d-csol)/csol*ecu
            eru = column.mat.eru
            if er>eru or er<-eru:    # control by reinforcement failure
                ecf = csol/(d-csol)*eru
            else:
                ecf = ecu
            phiu = ecf/csol
            Mu1 = N*(e0+l**2/np.pi**2*phiu)
            Mu2 = Msec
            return abs(Mu1-Mu2)
        x0candidate = np.linspace(0.01*Nccu, 0.99*Nccu, 50)
        dM = np.array([solvecolumn(Nui) for Nui in x0candidate])
        x0 = x0candidate[np.argmin(dM)]
        sol = minimize(solvecolumn, x0, bounds=((0.01*Nccu, 0.99*Nccu),))
        Ncol = sol.x[0]
        Mcol = fNM(Ncol)
        ci = fNci(Ncol)
        e1 = Mcol/Ncol - self.e0
        return Ncol, Mcol, e1


def getrccolumn(price):
    cost = Cost(price)
    # define material
    fc = 30.0
    fr = 300; Er = 200e3
    ecu = 0.0038
    eru = fr/Er*100
    Ec = 4700*np.sqrt(fc)
    fcp = 0.9*fc
    eco = 1.8*fcp/Ec
    def concss(ec):
        if ec<=eco:
            sc = fcp*(2*ec/eco-(ec/eco)**2)
        else:
            sc = fcp-0.15*fcp*(ec-eco)/(0.0038-eco)
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 700.
    Ac = np.pi/4*(h)**2
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([645., 1290, 1290, 1290, 1290, 1290, 645])
    xf = np.array([0.])
    xs = np.array([65., 103, 207.5, 350, 492.5, 597, 635])
    def bdist(x):
        if x<0 or x>h/2.:
            b = 0.
        else:
            b = 2.*np.sqrt((h/2.)**2-(x-(h/2.))**2)
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.axial_compression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.capacity(6000e3)
    print 'N = {} kN; M = {} kN-m; c={}'.format(N/1e3, M/1e6, csol)
    return column


def getfrpcolumn(price):
    cost = Cost(price)
    # define geometry
    h = 700.
    Ac = np.pi/4*h**2
    Afb = np.array([0.])
    tf = 6.0;
    Aft = np.pi*h*tf
    As = np.array([387, 774.,774.,774.,774.,774.,387.])
    xf = np.array([0.])
    xs = np.array([40., 81.5, 195, 350, 505, 618.5, 660])
    def bdist(x):
        if x<0 or x>h/2.:
            b = 0.
        else:
            b = 2.*np.sqrt((h/2.)**2-(x-(h/2.))**2)
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define material
    fc = 30.0
    fr = 500; Er = 45e3
    eco = 0.002
    psif=0.95; ke=0.55; ka=1.0; kb=1.0
    efe = ke*fr/Er
    fl = 2*Er*tf*efe/h
    print 'confinement ratio = {}'.format(fl/fc)
    fcc = fc+psif*3.3*ka*fl
    eccu = eco*(1.50+12*kb*fl/fc*(efe/eco)**0.45)
    print 'fcc = {}'.format(fcc)
    print 'ecc = {}'.format(eccu)
    if eccu>0.01: eccu=0.01
    eru = fr/Er*100
    Ec = 4700*np.sqrt(fc)
    def concss(ec):
        E2 = (fcc-fc)/eccu
        etp = 2*fc/(Ec-E2)
        if ec<=etp:
            sc = Ec*ec-(Ec-E2)**2/(4*fc)*ec**2
        elif ec>etp and ec<eccu:
            sc = fc+E2*ec
        else:
            sc = 0
        return sc
    ecu = eccu
    fc = fcc
    # reinss must be a vector function
    # reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    def reinss(x):
        sr = np.maximum(np.minimum(x*Er, fr), -fr)
        if isinstance(x,float) and sr<0:
            sr = 1e-6*sr
        else:
            sr[sr<0] = 1e-6*sr[sr<0]
        return sr
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.axial_compression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.capacity(7500.e3)
    print 'N = {} kN; M = {} kN-m; c={}'.format(N/1e3, M/1e6, csol)

    return column

if __name__ == '__main__':
    matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
            'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    price = Price(matprice=matprice)
    cost = Cost(price)
    ## RC example ex7.5 7.6 7.7 (Shu Shilin concrete textbook) checked
    # define material
    fc = 9.6
    fr = 300; Er = 200e3
    ecu = 0.0033
    eco = 0.002
    eru = fr/Er*100
    def concss(ec):    # must be scalar function
        if ec<=eco:
            sc = fc*(1.-(1.-ec/eco)**(2-20./60))
        else:
            sc = fc
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 500.
    Ac = 300.*500.
    Afb = np.array([0.])
    Aft = 0.
    As = np.array([1012., 308.])
    xf = np.array([0.])
    xs = np.array([40.,460.])
    def bdist(x):
        if x<0 or x>500:
            b = 0.
        else:
            b = 300.
        return b
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.axial_compression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.capacity(1150e3)
    print 'N = {} kN; M = {} kN-m'.format(N/1e3, M/1e6)

    ## RC example ex10.2 (Nelson Design of Reinforced Concrete)
    # define material
    fc = 27.6 #MPa: 4ksi
    fr = 413.7 #MPa: 60ksi
    Er = 200e3
    ecu = 0.003
    eco = (1-0.85)*ecu    # beta1 in Whitney block
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    # Whitney block
    def concss(ec, eco=eco, fc=fc):    # must be scalar function
        if ec<=eco:
            sc = 0.
        else:
            sc = 0.85*fc
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 609.6 #mm: 24 in
    b = 355.6 #mm: 14 in
    Ac = b*h
    Afb = np.array([0.])
    tft = 0.
    As = np.array([1935., 1935.])   # 3#9 for tension, 3#9 for compression
    xf = np.array([0.])
    xs = np.array([63.5,546.1])    # 2.5 in edge to rebar center
    def bdist(x, h=h, b=b):
        if x<0 or x>h:
            bsec = 0.
        else:
            bsec = b
        return bsec
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.axial_compression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.capacity(2774.36e3)    # 623.7k
    print 'N = {} kN; M = {} kN-m'.format(N/1e3, M/1e6)

    ## RC example ex10.6, a circular column (Nelson Design of Reinforced Concrete)
    # define material
    fc = 27.6 #MPa: 4ksi
    fr = 413.7 #MPa: 60ksi
    Er = 200e3
    ecu = 0.003
    eco = (1-0.85)*ecu    # beta1 in Whitney block
    # rupture strain of steel is assmued to be 100 times larger than yield strain
    eru = fr/Er*100
    # Whitney block
    def concss(ec, eco=eco, fc=fc):    # must be scalar function
        if ec<=eco:
            sc = 0.
        else:
            sc = 0.85*fc
        return sc
    # reinss must be a vector function
    reinss = lambda x: np.maximum(np.minimum(x*Er, fr), -fr)
    mat = Material(fc, fr, eco, ecu, eru, concss, reinss)
    # define geometry
    h = 508. #mm: 20 in
    Ac = np.pi*h**2/4.
    Afb = np.array([0.])
    tft = 0.
    Astotal = 0.004761281e6   # 8 rebars As (total) = 7.38 in2
    As = np.array([Astotal/8., Astotal/4., Astotal/4., Astotal/4., Astotal/8.])
    xf = np.array([0.])
    xs = np.array([63.5, 119.3, 254., 388.7, 444.5])    # 2.5in edge to rebar center
    def bdist(x, h=h, b=b):
        if x<0 or x>h:
            bsec = 0.
        else:
            bsec = 2*np.sqrt((h/2)**2-(x-h/2)**2)
        return bsec
    geo = Geometry('rc', h, Ac, Afb, Aft, As, xf, xs, bdist)
    # define beam
    column = Column(geo=geo, mat=mat, cost=cost)
    Nmax = column.axial_compression()
    print 'axial compression = {} kN'.format(Nmax/1e3)
    N,M,csol = column.capacity(3177.36e3)    # 714.3kip
    print 'N = {} kN; M = {} kN-m'.format(N/1e3, M/1e6)
    # textbook results: M=3428.64 kip-in = 387.384kNm
