import numpy as np

from section import Price

import sys
import os

discount = 0.02
mntcostdirect = {'Cwb': 11.5*9, 'Crm': 43.28*9, 'Cbc': 16.41*9, 'Ctp': 53.71*9}
mntcostindirect = {'Cruncar': 0.07*9, 'Cruntruck': 0.34*9, 'Cwage': 20.77*9, 'Cdriver': 24.54*9,
        'Ccargo': 3.64*9, 'Ttruck': 0.12, 'ADT': 8500., 'Ldetour': 2.9, 'Tdetour': 7, 'S': 50,
        'Ocars': 1.5, 'Otruck': 1.05}
mntplan = {'Mt': np.arange(20., 101., 20.), 'Me': 10.}
price = Price()

if __name__ == '__main__':
    test
