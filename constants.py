import numpy as np
from section import Price, Cost

# database parameters
HCOL = 1
L2HCOL = 2
TFTCOL = 3
EHCOL = 4
EHCOUPON = 6
FRCOL = 11
ERCOL = 12
RHOSCOL = 13
COVERCOL = 15
NLCOL = 17
FCCOL = 18
E0COL = 19
NUCOL = 20
MUCOL = 21
EHRUP = 23
FMCOL = 25

# design case parameter
HDESIGN = 400.0
LDESIGN = 1000.0
E0DESIGN = 0.15*HDESIGN
DS2DDESIGN = 0.80
FCODESIGN = 30.0
FRDESIGN = 400.0
ERDESIGN = 200e3
RHOSDESIGN = 0.01
NSDESIGN = 12.
EPHDESIGN = 1.5/100
EHDESIGN = 230e3
# parametric study
E2DARRAY = np.linspace(0.5, 0.25, 5)
RHOSARRAY = np.linspace(0.01, 0.05, 5)
DS2DARRAY = np.array([0.7, 0.8, 0.9])
FCC2FCOARRAY = np.array([1.25, 1.50, 1.75])
FCOARRAY = np.array([30., 40., 50.])


# code specific constants
ECO_GB = 0.002
ECU_GB = 0.0033

# rv parameters
H2MEAN = 1.524; HSTD = 6.35
FCO2MEAN = 1.20; FCOCOV = 0.10
FS2MEAN = 1.20; FSCOV = 0.10
FFRP2MEAN = 1.25; FFRPCOV = 0.12
EFRP2MEAN = 1.0; EFRPCOV = 0.10
L2MEAN = 1.00; LCOV = 0.25
D2MEAN = 1.05; DCOV = 0.10
# model error & correlation (subjected to change)
DNMEAN = 1.2479; DNSTD = 0.0989
DEMEAN = 1.0846; DESTD = 0.0389
rhodNE = 0.4434
rhoNM = 0.9580
# FRP efficiency factor (GB)
CFRP_K = 0.60

# partial safety factors
GAMMA_FC = 1.40
GAMMA_FS = 1.10
GAMMA_FF = 1.40
GAMMA_D_GB = 1.20
GAMMA_L_GB = 1.40
# strength reduciton factors
GAMMA_D_ACI = 1.40
GAMMA_L_ACI = 1.70

# irrelevant data for column initialization
matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
        'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
price = Price(matprice=matprice)
cost = Cost(price)
