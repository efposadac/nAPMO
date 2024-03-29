# file: timer.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from contextlib import contextmanager
import time
import operator
import napmo


class Timer(object):
    """
    Manages the execution times by context and for the global calculation
    """

    def __init__(self):
        super(Timer, self).__init__()
        self._start = time.time()
        self._end = 0.0
        self._blocks = {}

    @contextmanager
    def timeblock(self, label):
        """
        Calculates execution time for a context
        """
        self._blocks.setdefault(label, 0)
        start = time.time()
        try:
            yield
        finally:
            end = time.time()
            self._blocks[label] += (end - start) / napmo.threads

    def show_block(self, label):
        """
        Print information about certain context
        """
        print('\nTotal time for {0:s} : {1:5.4f} s'.format(
            label,
            self._blocks.get(label)))

    def show_summary(self):
        """
        Prints all context
        """
        print("""\n==================================================
Object:   {0:9s}
--------------------------------------------------
{1:<36s} {2:<12s}
--------------------------------------------------""".format(
            type(self).__name__,
            "Process",
            "Time (s)"))

        for item in sorted(self._blocks.items(), key=operator.itemgetter(1), reverse=True):
            print('{0:<36s} {1:<12.4f}'.format(
                item[0],
                item[1]))

        self._end = time.time()

        print("""--------------------------------------------------
{0:<36s} {1:<12.4f}
--------------------------------------------------
""".format(
            "Total time",
            (self._end - self._start) / napmo.threads))
