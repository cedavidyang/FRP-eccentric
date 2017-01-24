import numpy as np
import scipy.integrate as integrate

from section import Price, Cost, Material, Geometry
from beam import Beam, getrcbeam, getfrpbeam
from column import Column, getrccolumn, getfrpcolumn

import sys

class Bridge(object):
    def __init__(self, beam, beamno, column, colno, deck,
            span=None, width=None, height=None, distance=None):
        self.beam = beam
        self.column = column
        self.beamno = beamno
        self.colno = colno
        self.deck = deck
        self.span = span
        self.width = width
        self.height = height
        self.distance = distance

    def setcost(self, cost):
        self.cost = cost

    def setdistance(self, distance):
        seld.distance = distance

    def matcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        height = self.height
        beam = self.beam
        column = self.column
        cost = beamno*span*beam.matcost() + colno*height*column.matcost()
        return cost

    def transcost(self):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        height = self.height
        beam = self.beam
        column = self.column
        distance = self.distance
        cost = beamno*span*beam.transcost(distance) + colno*height*column.transcost(distance)
        return cost

    def mntcost(self, mntplan):
        beamno = self.beamno
        colno = self.colno
        span = self.span
        width = self.width
        height = self.height
        beam = self.beam
        column = self.column

        bplan = mntplan['beam']
        dplan = mntplan['deck']
        cplan = mntplan['column']

        # direct cost
        bcostpm = self.beam.directmntcost(bplan)
        dcostpm = self.deck.directmntcost(dplan)
        ccostpm = self.column.directmntcost(cplan)
        directcost = beamno*span*bcostpm + width*span*dcostpm + colno*height*ccostpm

        # indirect cost
        pindirect = self.beam.cost.price.mntprice['indirect']
        v = self.beam.cost.price.mntprice['discount']
        allmt = np.hstack((bplan['Mt'],dplan['Mt'],cplan['Mt']))
        allmt = np.unique(allmt)
        nmnt = allmt.size
        # price data
        Cruncar = pindirect['Cruncar']
        Cruntruck = pindirect['Cruntruck']
        Cwage = pindirect['Cwage']
        Cdriver = pindirect['Cdriver']
        Ccargo = pindirect['Ccargo']
        Ttruck = pindirect['Ttruck']
        Adt = pindirect['ADT']
        Ldetour = pindirect['Ldetour']
        Tdetour = pindirect['Tdetour']
        speed = pindirect['S']
        Ocar = pindirect['Ocar']
        Otruck = pindirect['Otruck']

        # indirect cost
        Crun = (Cruncar*(1-Ttruck)+Cruntruck*Ttruck)*Ldetour*Tdetour*Adt
        Ctime = (Cwage*Ocar*(1-Ttruck)+(Cdriver*Otruck+Ccargo)*Ttruck)*Ldetour*Tdetour*Adt/speed
        # Ctime = 0.
        indirectcost = (Crun+Ctime)
        # indirectcost = 0.
        # for t in allmt:
            # Crun = (Cruncar*(1-Ttruck)+Cruntruck*Ttruck)*Ldetour*Tdetour*Adt/(1+v)**t
            # Ctime = (Cwage*Ocar*(1-Ttruck)+(Cdriver*Otruck+Ccargo)*Ttruck)*Ldetour*Tdetour*Adt/speed/(1+v)**t
            # indirectcost += (Crun+Ctime)

        return (directcost, indirectcost)


if __name__ == '__main__':
    import copy
    # bridge distance
    distance = 10.
    # price information
    # matprice = {'Cconc': 945.-100, 'Csteel': 4855.74,
            # 'Cfb': 39810.23, 'Cft': 39810.23}
    # matprice = {'Cconc': 945.-100, 'Csteel': 3200.,
            # 'Cfb': 27812, 'Cft': 27812.}
    matprice = {'Cconc': 945.-100, 'Csteel': 5231.,
            'Cfb': 27812, 'Cft': 27812.}
    price = Price(matprice=matprice)
    rcbeam = getrcbeam(price)
    rcdeck = copy.deepcopy(rcbeam)
    rccolumn = copy.deepcopy(rcbeam)
    rcbridge = Bridge(rcbeam, 1, rccolumn, 0, rcdeck,
            span=8.0, width=0.0, height=6., distance=distance)

    # initial cost
    rcmatcost = rcbridge.matcost()
    rctranscost = rcbridge.transcost()
    print 'RC material cost = {}'.format(rcmatcost)
    print 'RC transportation cost = {}'.format(rctranscost)

    # maintenance cost
    rcbridge.beam.geo.addprop({'Arepair':400*1000.,'cover':40.})
    rcbridge.column.geo.addprop({'Arepair':np.pi*675*1000.,'cover':40.})
    Arepair = 1000.*1000.; cover=30.
    rcbridge.deck.geo.addprop({'Arepair':Arepair,'cover':cover})

    discount = 0.02
    mntpricedirect = {'Cwb': 11.5*8.72, 'Crm': 43.28*8.72, 'Cbc': 16.41*8.72, 'Ctp': 0.0*9}
    mntpriceindirect = {'Cruncar': 0.08*7.7, 'Cruntruck': 0.375*7.7, 'Cwage': 22.82*7.7, 'Cdriver': 26.97*7.7,
            'Ccargo': 4*7.7, 'Ttruck': 0.12, 'ADT': 8500., 'Ldetour': 2.9, 'Tdetour': 360., 'S': 50,
            'Ocar': 1.5, 'Otruck': 1.05}
    mntprice = {'direct':mntpricedirect, 'indirect':mntpriceindirect, 'discount':discount}

    rcbridge.beam.cost.price.setmntprice(mntprice)
    rcbridge.column.cost.price.setmntprice(mntprice)
    rcbridge.deck.cost.price.setmntprice(mntprice)

    # mntplan = {'Mt': np.arange(20., 101., 30.), 'Me': 20.}
    mntplan = {'Mt': np.arange(20., 101., 30.), 'Me': 20.}
    bridgemntplan = {'beam':mntplan, 'column':mntplan, 'deck':mntplan}
    mntplanN = {'Mt': None, 'Me': 10.}
    rcRarray = rcbridge.beam.lifecycleR(100, mntplan)
    rcRarrayN = rcbridge.beam.lifecycleR(100, mntplanN)
    rcmntcost = rcbridge.mntcost(bridgemntplan)
    print 'RC maintenace cost = {}'.format(rcmntcost)

    # FRP bridge
    # price information
    # matprice = {'Cconc': 884.3-100, 'Csteel': 4855.74,
            # 'Cfb': 39810.23, 'Cft': 39810.23}
    # matprice = {'Cconc': 884.3-100, 'Csteel': 3200.,
            # 'Cfb': 27812, 'Cft': 27812.}
    matprice = {'Cconc': 884.3-100, 'Csteel': 5231.,
            'Cfb': 28336., 'Cft': 28336.}
    price = Price(matprice=matprice)
    frpbeam = getfrpbeam(price)
    frpdeck = copy.deepcopy(frpbeam)
    frpcolumn = copy.deepcopy(frpbeam)
    frpbridge = Bridge(frpbeam, 1, frpcolumn, 0, frpdeck,
            span=8.0, width=0.0, height=6., distance=distance)
    # initial cost
    frpmatcost = frpbridge.matcost()
    frptranscost = frpbridge.transcost()
    print 'FRP material cost = {}'.format(frpmatcost)
    print 'FRP transportation cost = {}'.format(frptranscost)
    # maintenance cost
    mntplanN = {'Mt': None, 'Me': 10.}
    frpRarray = frpbeam.lifecycleR(100, mntplanN)
    import matplotlib.pyplot as plt
    plt.plot(np.arange(1,101), rcRarray)
    plt.plot(np.arange(1,101), rcRarrayN)
    plt.plot(np.arange(1,101), frpRarray)
    print "rc residual = {}".format(rcRarray[np.arange(1,101)==50])


    # # RC bridge pier
    # matprice = {'Cconc': 104.57*8.72, 'Csteel': 4871.,
            # 'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    # price = Price(matprice=matprice)
    # rccolumn = getrccolumn(price)
    # rcbridge = Bridge(rcbeam, 0, rccolumn, 1, rcdeck,
            # span=10., width=0.0, height=6., distance=1.)
    # # initial cost
    # rcmatcost = rcbridge.matcost()
    # rctranscost = rcbridge.transcost()
    # print 'RC material cost = {}'.format(rcmatcost)
    # print 'RC transportation cost = {}'.format(rctranscost)
    # # maintenance cost
    # rcbridge.column.geo.addprop({'Arepair':np.pi*675*1000.,'cover':40.})
    # rcbridge.column.cost.price.setmntprice(mntprice)
    # rcmntcost = rcbridge.mntcost(bridgemntplan)
    # print 'RC maintenace cost = {}'.format(rcmntcost)

    # # FRP bridge pier
    # matprice = {'Cconc': 104.57*8.72-100, 'Csteel': 4871.,
            # 'Cfb': 3.90e3*7.77, 'Cft': 3.90e3*7.77}
    # price = Price(matprice=matprice)
    # frpcolumn = getfrpcolumn(price)
    # # initial cost
    # frpmatcost = frpcolumn.matcost()
    # frptranscost = frpcolumn.transcost(1.)
    # print 'FRP material cost = {}'.format(frpmatcost)
    # print 'FRP transportation cost = {}'.format(frptranscost)
