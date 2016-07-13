# file: timer.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from contextlib import contextmanager
import time
import napmo


class Timer(object):
    """
    Manages the execution times by context and for the global calculation
    """

    def __init__(self):
        super(Timer, self).__init__()
        self._start = time.clock()
        self._end = 0.0
        self._blocks = {}

    @contextmanager
    def timeblock(self, label):
        self._blocks.setdefault(label, 0)
        start = time.clock()
        try:
            yield
        finally:
            end = time.clock()
            self._blocks[label] += (end - start) / napmo.threads

    def show_block(self, label):
        print('\nTotal time for {0:s} : {1:5.4f} s'.format(
            label,
            self._blocks.get(label)))

    def show_summary(self):
        print("""\n==================================================
Object:   {0:9s}
--------------------------------------------------
{1:<36s} {2:<12s}
--------------------------------------------------""".format(
            type(self).__name__,
            "Process",
            "Time (s)"))

        for label in self._blocks:
            print('{0:<36s} {1:<12.4f}'.format(
                label,
                self._blocks.get(label)))

        self._end = time.clock()

        print("""--------------------------------------------------
{0:<36s} {1:<12.4f}
--------------------------------------------------
""".format(
            "Total time",
            (self._end - self._start) / napmo.threads))
