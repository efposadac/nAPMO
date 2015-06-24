#!/usr/bin/env python3
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
import subprocess

from napmo.interfaces.stack import Stack


def gauss_chebishev_omp():
    nprocs = np.array([i for i in range(1, 9)])
    times = Stack()

    os.system('cd ../../src; make test')

    for nproc in nprocs:
        start_time = time.time()
        os.system(
                    'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../src; ' +
                    'export OMP_NUM_THREADS=' + str(nproc) +
                    '; ../../src/test')
        times.append(time.time() - start_time)

        print(nproc, times.peek())

    times = np.array(times)
    plt.title('Scaling in gauss_chebishev function')
    plt.plot(nprocs, times[0]/times, '-s', label='1e8 gauss_chebishev quad. 10 rep.')
    plt.xlabel('nprocs')
    plt.ylabel('speedup')
    plt.legend(loc=4)
    plt.savefig('gauss_chebishev.png')
    plt.close()


def mmi_omp():
    nprocs = np.array([i for i in range(1, 21)])
    times = Stack()

    for nproc in nprocs:
        # Calculate integral
        start_time = time.time()
        os.system(
                    'export OMP_NUM_THREADS=' + str(nproc) +
                    '; ./mmi_omp.py')
        times.append(time.time() - start_time)

        print(nproc, times.peek())

    times = np.array(times)
    plt.title('Scaling in MMI')
    plt.plot(nprocs, times[0]/times, '-s', label='C2 \int \\rho(r)')
    plt.xlabel('nprocs')
    plt.ylabel('speedup')
    plt.legend(loc=4)
    plt.savefig('mmi.png')
    plt.close()

if __name__ == '__main__':
    gauss_chebishev_omp()
    mmi_omp()
