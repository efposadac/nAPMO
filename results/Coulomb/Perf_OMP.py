#!/usr/bin/env python
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
import cProfile


def coulomb_omp():
    labels = ['1202 x 100', '1202 x 500', '1202 x 1000', '1202 x 2000']
    marks = ['-x', '-s', '-o', '-^']
    fig, ax = plt.subplots()
    for file in range(2, 6):
        nprocs = np.array([1, 2, 4, 8])
        times = []
        for nproc in nprocs:
            # Calculate integral
            start_time = time.time()
            os.system(
                'export OMP_NUM_THREADS=' + str(nproc) +
                '; ./Coulomb' + str(file) + '.py')
            times.append(time.time() - start_time)

            print(nproc, times[-1])

        times = np.array(times)

        plt.plot(nprocs, times[0] / times,
                 marks[file - 2], label=labels[file - 2])

        # plt.plot(nprocs, times,
        #          marks[file - 2], label=labels[file - 2])

        print(labels[file - 2], times)
    plt.xscale('log', basex=2)
    # plt.yscale('log')
    plt.xlabel('Number of Threads')
    plt.ylabel('Speedup')
    # plt.ylabel('Time (s)')
    plt.legend()
    ax.set_xticks(nprocs)
    ax.set_xticklabels(["1", "2", "4", "8"])
    plt.savefig('coulomb_speedup.png')
    plt.close()

if __name__ == '__main__':
    # gauss_chebishev_omp()
    coulomb_omp()
