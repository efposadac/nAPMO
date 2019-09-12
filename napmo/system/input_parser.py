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
import sys

def extract_keywork(data, beg=0):
    """
    Extract keywords from the input file starting in ``beg`` position

    Example::

        molecule {...} molecule is a keyword
        basis = ... basis is a keyword

    Args:
        data (str) : the input as string
        beg (int) : the starting position
    """
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
    """
    Converts a string in its corresponding data type.

    Example:

    .. doctest::

        >>> from napmo.system.input_parser import fix_casting
        >>> fix_casting('2')
        2
        >>> fix_casting('2.0')
        2.0
        >>> fix_casting('True')
        True
        >>> fix_casting('Test')
        'Test'
    """

    if not isinstance(item, str):
        return item

    try:
        return int(item)
    except ValueError:
        try:
            return np.float64(item)
        except ValueError:
            # Boolean
            if item == "True":
                return True
            if item == "False":
                return False
            # List
            if item.startswith('['):
                item = item.replace('[', '').replace(']', '').split(',')
                item = [fix_casting(val.strip()) for val in item]
                return item
            # Dictionary
            if item.find(':') >= 0:
                item = item.split(':')
                item = {item[0].strip(): fix_casting(item[1].strip())}
                return item
            return item


def raise_exception(e, error, message):
    try:
        raise e(error)
    except:
        print(message)
        raise


class InputParser(object):
    """
    Input parser class
    """

    def __init__(self, data):
        super(InputParser, self).__init__()

        action = {
            'molecule': self.load_molecule,
            'basis': self.load_basis,
            'scf': self.load_scf,
            'code': self.load_code,
            'grid': self.load_grid,
            'functional': self.load_functional
        }

        self.code = ''  # code to execute at the end of the calculation
        self.data = {}  # data of the molecule
        self.charges = {}  # options for charges and multiplicity
        self.scf = {}  # options for SCF engine
        self.grid = {}  # data for grids
        self.functional = {}  # data for functionals

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
                    raise_exception(
                        ValueError,
                        "Keywork action not known!",
                        'Check your input file ' + keyword + ' is undefined')

            pos += 1

    def load_molecule(self, data, group=True, options=None):
        """
        Load Molecule information from the input file (Only Cartesian coordinates for now)

        Args:
            data (str) : Relevant data from input corresponding to the keyword ``molecule``
            group (bool) : Whether the data is a group ie. ``{...}`` or just a variable
            options (str) : Additional options for the group ``[...]`` besides molecule
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
        """
        Loads the basis keyword

        Args:
            data (str) : Relevant data from input corresponding to the keyword ``basis``
            group (bool) : Whether the data is a group ie. ``{...}`` or just a variable ``basis = ...``

        """
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
        """
        Loads the scf keyword

        Args:
            data (str) : Relevant data from input corresponding to the keyword ``scf``
            group (bool) : Whether the data is a group ie. ``{...}`` or just a variable ``scf = ...``
        """
        options = {'uhf': {'method': 'uhf'},
                   'hf': {'method': 'rhf'},
                   'rhf': {'method': 'rhf'},
                   'analytic': {'kind': 'analytic'},
                   'numeric': {'kind': 'numeric'},
                   'direct': {'direct': True},
                   'tf': {'tf': True}
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

    def load_grid(self, data, group=True, options=None):

        aux = {}

        data = data.splitlines()
        
        for line in data:
            line = line.strip().split(' ')
            symbol = line[0].strip()

            if symbol not in self.data and symbol != 'e-':
                raise_exception(
                    ValueError,
                    "Particle not found!",
                    'Check the "grid" block in your input file ' + symbol + ' is undefined in "molecular" block')

            grid_spec = fix_casting(' '.join(line[1:]))

            if len(grid_spec) == 4:
                rtransform = napmo.PowerRadialTransform(
                    grid_spec[0], grid_spec[1], grid_spec[2])

                nrad = grid_spec[2]
                nang = grid_spec[3]

            elif len(grid_spec) == 3:
                rtransform = napmo.ChebyshevRadialTransform(
                    grid_spec[0], grid_spec[1])

                nrad = grid_spec[1]
                nang = grid_spec[2]

            else:
                rtransform = None
                nrad = grid_spec[0]
                nang = grid_spec[1]

            aux[symbol] = {'nrad': nrad,
                           'nang': nang,
                           'rtransform': rtransform}
            
        self.scf['grid'] = aux
        
    def load_functional(self, data, group=True, options=None):
        """
        Loads the functionals information for DFT calculations

        Args:
            data (str) : Relevant data from input corresponding to the keyword ``functional``
            group (bool) : Whether the data is a group ie. ``{...}`` or just a variable ``functional = ...``
        """

        aux = {}

        data = data.splitlines()
        
        for line in data:

            line = line.strip().split(' ')

            symbol = line[0].strip()
            if symbol not in self.data and symbol != 'e-':
                raise_exception(
                    ValueError,
                    "Particle not found!",
                    'Check the "functional" block in your input file ' + symbol + ' is undefined in "molecular" block')

            if len(line) > 3:
                raise_exception(
                    ValueError,
                    "There is an extra term in a functional line!",
                    'Check the "functional" block in your input file')

            elif len(line) == 3:
                otherSymbol= line[1].strip()

                if otherSymbol not in self.data and symbol != 'e-':
                    raise_exception(
                    ValueError, 
                    "Particle not found!",
                    'Check the "functional" block in your input file ' + symbol + ' is undefined in "molecular" block')

                functional= line[2].strip()
                aux[symbol+'//'+otherSymbol]= functional

            else:
                functional= line[1].strip()
                aux[symbol] = functional
                
        self.functional = aux
        print self.functional

        
    def load_code(self, data, group=True, options=None):
        """
        Loads the code keyword

        Args:
            data (str) : Relevant data from input corresponding to the keyword ``code``
        """
        self.code = data
