# file: Perf_OMP.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import matplotlib.pyplot as plt
import numpy as np
import time
import os
import sys


def gauss_chebishev_omp():
    nprocs = np.array([i for i in range(1, 9)])
    times = []

    for nproc in nprocs:
        start_time = time.time()
        os.system(
                    'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../src; ' +
                    'export OMP_NUM_THREADS=' + str(nproc) +
                    '; ../../src/test')
        times.append(time.time() - start_time)

    times = np.array(times)
    print(times)
    plt.title('Scaling in gauss_chebishev function')
    plt.plot(nprocs, times[0]/times, '-s', label='1e8 gauss_chebishev quadrature calculation x 10 times')
    plt.xlabel('nprocs')
    plt.ylabel('time (s)')
    plt.legend()
    plt.savefig('gauss_chebishev.png')
    plt.close()

if __name__ == '__main__':
    gauss_chebishev_omp()
