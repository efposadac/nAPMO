# file: input_parser.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import napmo
import numpy as np
import os
import re


def extract_keywork(data, beg=0):

    keyword = None
    group = False
    ingroup = False
    options = None

    pos = min(data.find('=', beg), data.find(
        '{', beg), key=lambda x: 1e100 if x == -1 else x)

    eof = data.find('{', pos) == -1
    sof = data.rfind('\n', 0, pos) == -1

    if pos >= 0:
        if sof:
            keyword = data[0:pos].strip(' \n')
        else:
            keyword = data[data.rfind('\n', 0, pos):pos].strip(' \n')

        if keyword:
            if keyword.find('[') > -1:
                options = keyword[keyword.find(
                    '['):keyword.find(']')].strip(' []')
                keyword = keyword[:keyword.find('[')].strip()

        group = data.find('{', pos, data.find('\n', pos)) >= 0

        if eof:
            ingroup = not group and (data.find('{', pos) < data.find('}', pos))
        else:
            ingroup = not group and (data.find('{', pos) > data.find('}', pos))

    return keyword, pos, options, group, ingroup


def fix_casting(item):
    try:
        return int(item)
    except ValueError:
        try:
            return float(item)
        except ValueError:
            if item == "True":
                return True
            if item == "False":
                return False
            return item


def raise_exception(e, error, message):
    try:
        raise e(error)
    except:
        print(message)
        raise


class InputParser(object):
    """
    Input parser
    """

    def __init__(self, data):
        super(InputParser, self).__init__()

        action = {
            'molecule': self.load_molecule,
            'basis': self.load_basis,
            'scf': self.load_scf
        }

        self.var = []  # variables for some math in input
        self.data = {}  # data of the molecule
        self.charges = {}  # options for charges and multiplicity
        self.scf = {}  # options for SCF engine

        # remove all comments from data and blank lines
        data = re.sub('""".*\n?', '', data)
        data = re.sub('#.*\n?', '', data)
        data = re.sub('(?imu)^\s*\n', '', data)

        pos = 0
        while True:
            keyword, pos, options, group, ingroup = extract_keywork(data, pos)

            if pos < 0 or not keyword:
                break

            if group:
                keydata = data[
                    data.find('{', pos) + 1: data.find('}', pos)].strip(' \n')

                action[keyword](keydata, group, options)

            if not group and not ingroup:

                if keyword in action:
                    keydata = data[pos:data.find('\n', pos)].strip('= \n')

                    action[keyword](keydata, ingroup)
                else:
                    self.var.append(
                        keyword + data[pos:data.find('\n', pos)].strip()
                    )

            pos += 1

    def load_molecule(self, data, group=True, options=None):
        """
        Load Molecule information from the input file (Only Cartesian coordinates for now)
        """
        atomic_data = napmo.AtomicElementsDatabase()
        keys = ['charge', 'multiplicity']

        data = data.splitlines()
        for line in data:
            line = ' '.join(line.split())
            line = line.strip().split(' ')

            symbol = line[0].strip()
            self.data.setdefault(symbol, [])

            origin = np.array(line[1:], dtype=np.float64)

            particle = {}
            particle.setdefault('origin', origin)
            particle.setdefault('basis_name', None)
            particle.setdefault('basis_file', None)

            if symbol in atomic_data:
                particle.setdefault('quantum', False)
            else:
                particle.setdefault('quantum', True)

            self.data[symbol].append(particle)

            aux = symbol.find('_')
            if aux > 0:  # add atom
                symbol = symbol[:aux]
                self.data.setdefault(symbol, [])

                atom = {}
                atom.setdefault('origin', origin)
                atom.setdefault('basis_name', None)
                atom.setdefault('basis_file', None)
                atom.setdefault('quantum', True)

                self.data[symbol].append(atom)

        if options:
            options = ' '.join(options.split())
            self.charges = {aux.split(':')[0].strip(): {
                keys[i - 1]: fix_casting(aux.split(':')[i].strip())
                for i in range(1, len(aux.split(':')))}
                for aux in options.split(',')}

    def load_basis(self, data, group=True, options=None):
        atomic_data = napmo.AtomicElementsDatabase()
        data = data.splitlines()

        if len(self.data) < 1:
            raise_exception(
                ValueError,
                "Molecular block undefined",
                'Modify your input file and define "molecular" before "basis"')

        if group:
            for line in data:
                line = ' '.join(line.split())
                line = line.strip().split(' ')

                symbol = line[0].strip()

                if symbol not in self.data and symbol != 'e-':
                    raise_exception(
                        ValueError,
                        "Particle not found!",
                        'Check the "basis" block in your input file ' + symbol + ' is undefined in "molecular" block')

                basis_name = line[1].strip()
                basis_file = os.path.join(napmo.basis_dir, basis_name)

                if len(line) > 2:
                    basis_file = line[2].strip()

                if symbol == 'e-':
                    for key in self.data:
                        if key in atomic_data:
                            for particle in self.data[key]:
                                particle['basis_name'] = basis_name
                                particle['basis_file'] = basis_file
                else:
                    for particle in self.data[symbol]:
                        particle['basis_name'] = basis_name
                        particle['basis_file'] = basis_file

        else:
            species = set()
            for key in self.data:
                if key in atomic_data:
                    species.add('e-')
                    continue
                for particle in self.data[key]:
                    if particle.get('quantum'):
                        species.add(key)

            if len(species) > 1:
                raise_exception(
                    ValueError,
                    "More than one species found!",
                    'Modify your input file and define a "basis" for each quantum species')

            basis_name = data[-1].strip()
            basis_file = os.path.join(napmo.basis_dir, basis_name)

            for key in self.data:
                for particle in self.data[key]:
                    particle['basis_name'] = basis_name
                    particle['basis_file'] = basis_file

    def load_scf(self, data, group=True, options=None):
        options = {'uhf': {'method': 'uhf'},
                   'hf': {'method': 'rhf'},
                   'rhf': {'method': 'rhf'},
                   'analytic': {'kind': 'analytic'},
                   'numeric': {'kind': 'numeric'},
                   'direct': {'direct': True}
                   }

        aux = {}
        data = data.splitlines()
        for line in data:
            line = ' '.join(line.split())
            line = line.strip()
            keyword, pos = extract_keywork(line)[0:2]

            if keyword:
                value = line[pos + 1:].strip()
                self.scf.setdefault(keyword, fix_casting(value))
            else:
                aux.update(options[line])

        self.scf.update(aux)
