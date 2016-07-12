import re

lvalue = {'0': 's', '1': 'p', '2': 'd', '3': 'f', '4': 'g'}

symbol = 'H'
name = 'hydrogen'
basis_data = []

aux = {}

with open('STO-3G', 'r') as f:
    data = f.read()

data = re.sub('#.*\n?', '', data)
data = re.sub('(?imu)^\s*\n', '', data)

data = data[data.find('O-' + name.upper()):].splitlines()

line_index = 2
for cont in range(int(data[1].strip())):
    line = data[line_index].strip().split()

    aux['angular'] = lvalue.get(line[1].strip())
    aux['cont'] = []
    aux['prim'] = []

    line_index += 1

    for prim in range(int(line[2].strip())):
        line = data[line_index].strip().split()

        aux['prim'].append(float(line[0].strip()))
        aux['cont'].append(float(line[1].strip()))

    basis_data.append(aux)

print(basis_data)
