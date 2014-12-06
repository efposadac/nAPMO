# file: databases.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

AtomicElementsDatabase = {

'H':

{'name':'Hydrogen','symbol':'H','atomicNumber':1,'mass':1.00794,
'melting_point':13.81,'boiling_point':20.28,'density':0.084,
'electron_affinity':-73,'ionization_energy_1':1312,'electronegativity':2.1,
'covalent_radius':0.30,'atomic_radius':0.25,'vanderwaals_radius':1.20},

'isotopes_H':{
1:{'mass_number':1,'atomicWeight':1.0078250321,'abundance':99.9885,'nuclearSpin':0.5},
2:{'mass_number':2,'atomicWeight':2.0141017780,'abundance':0.0115,'nuclearSpin':1},
3:{'mass_number':3,'atomicWeight':3.0160492675,'abundance':0.0,'nuclearSpin':0.5},
4:{'mass_number':4,'atomicWeight':4.02783,'abundance':0,'nuclearSpin':-2.0},
5:{'mass_number':5,'atomicWeight':5.03954,'abundance':0,'nuclearSpin':0},
6:{'mass_number':6,'atomicWeight':6.04494,'abundance':0,'nuclearSpin':0},
'most_abundant':1},

'He':

{'name':'Helium','symbol':'He','atomicNumber':2,'mass':4.002602,
'melting_point':0.95,'boiling_point':4.216,'density':0.1785,
'electron_affinity':21,'ionization_energy_1':2375.1,
'covalent_radius':1.22,'atomic_radius':31,'vanderwaals_radius':140,'electronegativity':0},

'isotopes_He':{
3:{'mass_number':3,'atomicWeight':3.0160293097,'abundance':0.000137,'nuclearSpin':0.5},
4:{'mass_number':4,'atomicWeight':4.0026032497,'abundance':99.999863,'nuclearSpin':0.0},
5:{'mass_number':5,'atomicWeight':5.012220,'abundance':0,'nuclearSpin':-1.5},
6:{'mass_number':6,'atomicWeight':6.0188881,'abundance':0,'nuclearSpin':0.0},
7:{'mass_number':7,'atomicWeight':7.028030,'abundance':0,'nuclearSpin':-1.5},
8:{'mass_number':8,'atomicWeight':8.033922,'abundance':0,'nuclearSpin':0.0},
9:{'mass_number':9,'atomicWeight':9.043820,'abundance':0,'nuclearSpin':-0.5},
10:{'mass_number':10,'atomicWeight':10.052400,'abundance':0,'nuclearSpin':0.0},
'most_abundant':4},

'Li':

{'name':'Lithium','symbol':'Li','atomicNumber':3,'mass':6.941,
'melting_point':453.7,'boiling_point':1615,'density':0.53,
'electron_affinity':-60,'ionization_energy_1':520.867,
'covalent_radius':1.23,'atomic_radius':1.45,'vanderwaals_radius':1.82,'electronegativity':0.98},

'isotopes_Li':{
4:{'mass_number':4,'atomicWeight':4.02718,'abundance':0,'nuclearSpin':-2.0},
5:{'mass_number':5,'atomicWeight':5.012540,'abundance':0,'nuclearSpin':-1.5},
6:{'mass_number':6,'atomicWeight':6.0151223,'abundance':7.59,'nuclearSpin':1.0},
7:{'mass_number':7,'atomicWeight':7.0160040,'abundance':92.41,'nuclearSpin':-1.5},
8:{'mass_number':8,'atomicWeight':8.0224867,'abundance':0,'nuclearSpin':2.0},
9:{'mass_number':9,'atomicWeight':9.0267891,'abundance':0,'nuclearSpin':-1.5},
10:{'mass_number':10,'atomicWeight':10.035481,'abundance':0,'nuclearSpin':0},
11:{'mass_number':11,'atomicWeight':11.043796,'abundance':0,'nuclearSpin':-1.5},
12:{'mass_number':12,'atomicWeight':12.05378,'abundance':0,'nuclearSpin':0},
'most_abundant':7},

'Be':

{'name':'Beryllium','symbol':'Be','atomicNumber':4,'standarAtomicWeight':9.012182,
'melting_point':453.7,'boiling_point':1615,'density':0.53,
'electron_affinity':-60,'ionization_energy_1':520.867,
'covalent_radius':0.89,'atomic_radius':1.45,'vanderwaals_radius':1.82,'electronegativity':0.98},

'isotopes_Be':{
5:{'mass_number':5,'atomicWeight':5.04079,'abundance':0,'nuclearSpin':0},
6:{'mass_number':6,'atomicWeight':6.019726,'abundance':0,'nuclearSpin':0},
7:{'mass_number':7,'atomicWeight':7.0169292,'abundance':0,'nuclearSpin':0},
8:{'mass_number':8,'atomicWeight':8.00530509,'abundance':0,'nuclearSpin':0},
9:{'mass_number':9,'atomicWeight':9.0121821,'abundance':100,'nuclearSpin':0},
10:{'mass_number':10,'atomicWeight':10.0135337,'abundance':0,'nuclearSpin':0},
11:{'mass_number':11,'atomicWeight':11.021658,'abundance':0,'nuclearSpin':0},
12:{'mass_number':12,'atomicWeight':12.026921,'abundance':0,'nuclearSpin':0},
13:{'mass_number':13,'atomicWeight':13.03613,'abundance':0,'nuclearSpin':0},
14:{'mass_number':14,'atomicWeight':14.04282,'abundance':0,'nuclearSpin':0},
'most_abundant':9},

'B':

{'name':'Boron','symbol':'B','atomicNumber':5,'standarAtomicWeight':10.811,
'melting_point':2365,'boiling_point':4275,'density':2.46,
'electron_affinity':-27,'ionization_energy_1':801.587,
'covalent_radius':0.88,'atomic_radius':85,'vanderwaals_radius':0,'electronegativity':2.04},

'isotopes_B':{
7:{'mass_number':7,'atomicWeight':7.029920,'abundance':0,'nuclearSpin':0},
8:{'mass_number':8,'atomicWeight':8.0246067,'abundance':0,'nuclearSpin':0},
9:{'mass_number':9,'atomicWeight':9.0133288,'abundance':0,'nuclearSpin':0},
10:{'mass_number':10,'atomicWeight':10.0129370,'abundance':19.9,'nuclearSpin':0},
11:{'mass_number':11,'atomicWeight':11.0093055,'abundance':80.1,'nuclearSpin':0},
12:{'mass_number':12,'atomicWeight':12.0143521,'abundance':0,'nuclearSpin':0},
13:{'mass_number':13,'atomicWeight':13.0177803,'abundance':0,'nuclearSpin':0},
14:{'mass_number':14,'atomicWeight':14.025404,'abundance':0,'nuclearSpin':0},
15:{'mass_number':15,'atomicWeight':15.031097,'abundance':0,'nuclearSpin':0},
16:{'mass_number':16,'atomicWeight':16.039810,'abundance':0,'nuclearSpin':0},
17:{'mass_number':17,'atomicWeight':17.04693,'abundance':0,'nuclearSpin':0},
18:{'mass_number':18,'atomicWeight':18.05617,'abundance':0,'nuclearSpin':0},
19:{'mass_number':19,'atomicWeight':19.06373,'abundance':0,'nuclearSpin':0},
'most_abundant':11},

'C':

{'name':'Carbon','symbol':'C','atomicNumber':6,'mass':12.0107,
'melting_point':3825,'boiling_point':5100,'density':3.51,
'electron_affinity':-122,'ionization_energy_1':1087.72,
'covalent_radius':0.77,'atomic_radius':0.70,'vanderwaals_radius':1.70,'electronegativity':2.55},

'isotopes_C':{
8:{'mass_number':8,'atomicWeight':8.037675,'abundance':0,'nuclearSpin':0},
9:{'mass_number':9,'atomicWeight':9.0310401,'abundance':0,'nuclearSpin':0},
10:{'mass_number':10,'atomicWeight':10.0168531,'abundance':0,'nuclearSpin':0},
11:{'mass_number':11,'atomicWeight':11.0114338,'abundance':0,'nuclearSpin':0},
12:{'mass_number':12,'atomicWeight':12.0000000,'abundance':98.93,'nuclearSpin':0},
13:{'mass_number':13,'atomicWeight':13.0033548378,'abundance':1.07,'nuclearSpin':0},
14:{'mass_number':14,'atomicWeight':14.003241988,'abundance':0,'nuclearSpin':0},
15:{'mass_number':15,'atomicWeight':15.0105993,'abundance':0,'nuclearSpin':0},
16:{'mass_number':16,'atomicWeight':16.014701,'abundance':0,'nuclearSpin':0},
17:{'mass_number':17,'atomicWeight':17.022584,'abundance':0,'nuclearSpin':0},
18:{'mass_number':18,'atomicWeight':18.026760,'abundance':0,'nuclearSpin':0},
19:{'mass_number':19,'atomicWeight':19.03525,'abundance':0,'nuclearSpin':0},
20:{'mass_number':20,'atomicWeight':20.04032,'abundance':0,'nuclearSpin':0},
21:{'mass_number':21,'atomicWeight':21.04934,'abundance':0,'nuclearSpin':0},
22:{'mass_number':22,'atomicWeight':22.05645,'abundance':0,'nuclearSpin':0},
'most_abundant':12},

'N':

{'name':'Nitrogen','symbol':'N','atomicNumber':7,'mass':14.0067,
'melting_point':63.15,'boiling_point':77.344,'density':1.17,
'electron_affinity':7,'ionization_energy_1':1403.98,
'covalent_radius':0.70,'atomic_radius':0.65,'vanderwaals_radius':1.55,'electronegativity':3.04},

'isotopes_N':{
10:{'mass_number':10,'atomicWeight':10.04262,'abundance':0,'nuclearSpin':0},
11:{'mass_number':11,'atomicWeight':11.02680,'abundance':0,'nuclearSpin':0},
12:{'mass_number':12,'atomicWeight':12.0186132,'abundance':0,'nuclearSpin':1.0},
13:{'mass_number':13,'atomicWeight':13.00573858,'abundance':0,'nuclearSpin':-0.5},
14:{'mass_number':14,'atomicWeight':14.0030740052,'abundance':99.632,'nuclearSpin':1.0},
15:{'mass_number':15,'atomicWeight':15.0001088984,'abundance':0.368,'nuclearSpin':-0.5},
16:{'mass_number':16,'atomicWeight':16.0061014,'abundance':0,'nuclearSpin':-2.0},
17:{'mass_number':17,'atomicWeight':17.008450,'abundance':0,'nuclearSpin':-0.5},
18:{'mass_number':18,'atomicWeight':18.014082,'abundance':0,'nuclearSpin':-1},
19:{'mass_number':19,'atomicWeight':19.017027,'abundance':0,'nuclearSpin':0},
20:{'mass_number':20,'atomicWeight':20.023370,'abundance':0,'nuclearSpin':0},
21:{'mass_number':21,'atomicWeight':21.02709,'abundance':0,'nuclearSpin':0},
22:{'mass_number':22,'atomicWeight':22.03444,'abundance':0,'nuclearSpin':0},
23:{'mass_number':23,'atomicWeight':23.04051,'abundance':0,'nuclearSpin':0},
24:{'mass_number':24,'atomicWeight':24.05050,'abundance':0,'nuclearSpin':0},
'most_abundant':14},

'O':

{'name':'Oxygen','symbol':'O','atomicNumber':8,'mass':15.9994,
'melting_point':54.8,'boiling_point':90.188,'density':1.33,
'electron_affinity':-141,'ionization_energy_1':1315.5,
'covalent_radius':0.66,'atomic_radius':0.60,'vanderwaals_radius':1.52,'electronegativity':3.44},

'isotopes_O':{
12:{'mass_number':12,'atomicWeight':12.034405,'abundance':0,'nuclearSpin':0},
13:{'mass_number':13,'atomicWeight':13.024810,'abundance':0,'nuclearSpin':0},
14:{'mass_number':14,'atomicWeight':14.00859529,'abundance':0,'nuclearSpin':0},
15:{'mass_number':15,'atomicWeight':15.0030654,'abundance':0,'nuclearSpin':0},
16:{'mass_number':16,'atomicWeight':15.9949146221,'abundance':99.757,'nuclearSpin':0},
17:{'mass_number':17,'atomicWeight':16.99913150,'abundance':0.038,'nuclearSpin':0},
18:{'mass_number':18,'atomicWeight':17.9991604,'abundance':0.205,'nuclearSpin':0},
19:{'mass_number':19,'atomicWeight':19.003579,'abundance':0,'nuclearSpin':0},
20:{'mass_number':20,'atomicWeight':20.0040762,'abundance':0,'nuclearSpin':0},
21:{'mass_number':21,'atomicWeight':21.008655,'abundance':0,'nuclearSpin':0},
22:{'mass_number':22,'atomicWeight':22.009970,'abundance':0,'nuclearSpin':0},
23:{'mass_number':23,'atomicWeight':23.01569,'abundance':0,'nuclearSpin':0},
24:{'mass_number':24,'atomicWeight':24.02037,'abundance':0,'nuclearSpin':0},
25:{'mass_number':25,'atomicWeight':25.02914,'abundance':0,'nuclearSpin':0},
26:{'mass_number':26,'atomicWeight':26.03775,'abundance':0,'nuclearSpin':0},
'most_abundant':16},

'F':

{'name':'Fluorine','symbol':'F','atomicNumber':9,'mass':18.9984032,
'melting_point':53.55,'boiling_point':85,'density':1.58,
'electron_affinity':-328,'ionization_energy_1':1682.97,
'covalent_radius':0.58,'atomic_radius':50,'vanderwaals_radius':147,'electronegativity':3.98},

'isotopes_F':{
14:{'mass_number':14,'atomicWeight':14.03608,'abundance':0,'nuclearSpin':0},
15:{'mass_number':15,'atomicWeight':15.01801,'abundance':0,'nuclearSpin':0},
16:{'mass_number':16,'atomicWeight':16.011466,'abundance':0,'nuclearSpin':0},
17:{'mass_number':17,'atomicWeight':17.00209524,'abundance':0,'nuclearSpin':0},
18:{'mass_number':18,'atomicWeight':18.0009377,'abundance':0,'nuclearSpin':0},
19:{'mass_number':19,'atomicWeight':18.99840320,'abundance':100,'nuclearSpin':0},
20:{'mass_number':20,'atomicWeight':19.99998132,'abundance':0,'nuclearSpin':0},
21:{'mass_number':21,'atomicWeight':20.9999489,'abundance':0,'nuclearSpin':0},
22:{'mass_number':22,'atomicWeight':22.002999,'abundance':0,'nuclearSpin':0},
23:{'mass_number':23,'atomicWeight':23.003570,'abundance':0,'nuclearSpin':0},
24:{'mass_number':24,'atomicWeight':24.008100,'abundance':0,'nuclearSpin':0},
25:{'mass_number':25,'atomicWeight':25.012090,'abundance':0,'nuclearSpin':0},
26:{'mass_number':26,'atomicWeight':26.01963,'abundance':0,'nuclearSpin':0},
27:{'mass_number':27,'atomicWeight':27.02689,'abundance':0,'nuclearSpin':0},
28:{'mass_number':28,'atomicWeight':28.03567,'abundance':0,'nuclearSpin':0},
29:{'mass_number':29,'atomicWeight':29.04326,'abundance':0,'nuclearSpin':0},
'most_abundant':19},

'Ne':

{'name':'Neon','symbol':'Ne','atomicNumber':10,'mass':20.1797,
'melting_point':24.55,'boiling_point':27.1,'density':0.8999,
'electron_affinity':29,'ionization_energy_1':2083.08,
'covalent_radius':1.60,'atomic_radius':38,'vanderwaals_radius':154,'electronegativity':0},

'isotopes_Ne':{
16:{'mass_number':16,'atomicWeight':16.025757,'abundance':0,'nuclearSpin':0},
17:{'mass_number':17,'atomicWeight':17.017700,'abundance':0,'nuclearSpin':0},
18:{'mass_number':18,'atomicWeight':18.0056971,'abundance':0,'nuclearSpin':0},
19:{'mass_number':19,'atomicWeight':19.0018798,'abundance':0,'nuclearSpin':0},
20:{'mass_number':20,'atomicWeight':19.9924401759,'abundance':90.48,'nuclearSpin':0},
21:{'mass_number':21,'atomicWeight':20.99384674,'abundance':0.27,'nuclearSpin':0},
22:{'mass_number':22,'atomicWeight':21.99138551,'abundance':9.25,'nuclearSpin':0},
23:{'mass_number':23,'atomicWeight':22.99446734,'abundance':0,'nuclearSpin':0},
24:{'mass_number':24,'atomicWeight':23.993615,'abundance':0,'nuclearSpin':0},
25:{'mass_number':25,'atomicWeight':24.997790,'abundance':0,'nuclearSpin':0},
26:{'mass_number':26,'atomicWeight':26.000460,'abundance':0,'nuclearSpin':0},
27:{'mass_number':27,'atomicWeight':27.00762,'abundance':0,'nuclearSpin':0},
28:{'mass_number':28,'atomicWeight':28.01211,'abundance':0,'nuclearSpin':0},
29:{'mass_number':29,'atomicWeight':29.01935,'abundance':0,'nuclearSpin':0},
30:{'mass_number':30,'atomicWeight':30.02387,'abundance':0,'nuclearSpin':0},
31:{'mass_number':31,'atomicWeight':31.03311,'abundance':0,'nuclearSpin':0},
32:{'mass_number':32,'atomicWeight':32.03991,'abundance':0,'nuclearSpin':0},
'most_abundant':20},

'Na':

{'name':'Sodium','symbol':'Na','atomicNumber':11,'mass':22.98977,
'melting_point':371.0,'boiling_point':1156.0,'density':0.97,
'electron_affinity':-53.0,'ionization_energy_1':496.4274,
'covalent_radius':1.66,'atomic_radius':1.80,'vanderwaals_radius':2.27,'electronegativity':0.93},

'isotopes_Na':{
23:{'mass_number':23,'atomicWeight':22.98976967,'abundance':100,'nuclearSpin':0},
'most_abundant':23},

'Mg':

{'name':'Magnesium','symbol':'Mg','atomicNumber':12,'mass':24.3050,
'melting_point':922,'boiling_point':1380,'density':1.74,
'electron_affinity':19,'ionization_energy_1':738.6036,
'covalent_radius':1.36,'atomic_radius':1.50,'vanderwaals_radius':1.73,'electronegativity':1.31},

'isotopes_Mg':{
24:{'mass_number':24,'atomicWeight':23.98504190,'abundance':78.99,'nuclearSpin':0},
25:{'mass_number':25,'atomicWeight':24.98583702,'abundance':10.00,'nuclearSpin':0},
26:{'mass_number':26,'atomicWeight':25.98259304,'abundance':11.01,'nuclearSpin':0},
'most_abundant':24},

'Al':

{'name':'Aluminum','symbol':'Al','atomicNumber':13,'mass':26.981538,
'melting_point':933.5,'boiling_point':2740,'density':2.7,
'electron_affinity':-43.0,'ionization_energy_1':578.2476,
'covalent_radius':1.25,'atomic_radius':1.25,'vanderwaals_radius':0,'electronegativity':1.61},

'isotopes_Al':{
27:{'mass_number':27,'atomicWeight':26.98153844,'abundance':100,'nuclearSpin':0},
'most_abundant':27},

'Si':

{'name':'Silicon','symbol':'Si','atomicNumber':14,'mass':28.0855,
'melting_point':1683,'boiling_point':2630,'density':2.33,
'electron_affinity':-134.0,'ionization_energy_1':787.3866,
'covalent_radius':1.17,'atomic_radius':1.10,'vanderwaals_radius':2.10,'electronegativity':1.91},

'isotopes_Si':{
28:{'mass_number':28,'atomicWeight':27.9769265327,'abundance':92.2297,'nuclearSpin':0},
29:{'mass_number':29,'atomicWeight':28.97649472,'abundance':4.6832,'nuclearSpin':0},
30:{'mass_number':30,'atomicWeight':29.97377022,'abundance':3.0872,'nuclearSpin':0},
'most_abundant':28},

'P':

{'name':'Phosphorus','symbol':'P','atomicNumber':15,'mass':30.973761,
'melting_point':317.3,'boiling_point':553.0,'density':1.82,
'electron_affinity':-72.0,'ionization_energy_1':1012.7496,
'covalent_radius':1.10,'atomic_radius':1.00,'vanderwaals_radius':1.80,'electronegativity':2.19},

'isotopes_P':{
31:{'mass_number':31,'atomicWeight':30.97376151,'abundance':100,'nuclearSpin':0},
'most_abundant':31},

'S':

{'name':'Sulfur','symbol':'S','atomicNumber':16,'mass':32.065,
'melting_point':392.2,'boiling_point':717.82,'density':2.06,
'electron_affinity':-200,'ionization_energy_1':1000.78,
'covalent_radius':1.04,'atomic_radius':1.00,'vanderwaals_radius':1.80,'electronegativity':2.58},

'isotopes_S':{
26:{'mass_number':26,'atomicWeight':26.02788,'abundance':0,'nuclearSpin':0},
27:{'mass_number':27,'atomicWeight':27.01880,'abundance':0,'nuclearSpin':0},
28:{'mass_number':28,'atomicWeight':28.00437,'abundance':0,'nuclearSpin':0},
29:{'mass_number':29,'atomicWeight':28.996610,'abundance':0,'nuclearSpin':0},
30:{'mass_number':30,'atomicWeight':29.984903,'abundance':0,'nuclearSpin':0},
31:{'mass_number':31,'atomicWeight':30.9795544,'abundance':0,'nuclearSpin':0},
32:{'mass_number':32,'atomicWeight':31.97207069,'abundance':94.93,'nuclearSpin':0},
33:{'mass_number':33,'atomicWeight':32.97145850,'abundance':0.76,'nuclearSpin':0},
34:{'mass_number':34,'atomicWeight':33.96786683,'abundance':4.29,'nuclearSpin':0},
35:{'mass_number':35,'atomicWeight':34.96903214,'abundance':0,'nuclearSpin':0},
36:{'mass_number':36,'atomicWeight':35.96708088,'abundance':0.02,'nuclearSpin':0},
37:{'mass_number':37,'atomicWeight':36.97112572,'abundance':0,'nuclearSpin':0},
38:{'mass_number':38,'atomicWeight':37.971163,'abundance':0,'nuclearSpin':0},
39:{'mass_number':39,'atomicWeight':38.975140,'abundance':0,'nuclearSpin':0},
40:{'mass_number':40,'atomicWeight':39.97547,'abundance':0,'nuclearSpin':0},
41:{'mass_number':41,'atomicWeight':40.98003,'abundance':0,'nuclearSpin':0},
42:{'mass_number':42,'atomicWeight':41.98149,'abundance':0,'nuclearSpin':0},
43:{'mass_number':43,'atomicWeight':42.98660,'abundance':0,'nuclearSpin':0},
44:{'mass_number':44,'atomicWeight':43.98832,'abundance':0,'nuclearSpin':0},
45:{'mass_number':45,'atomicWeight':44.99482,'abundance':0,'nuclearSpin':0},
46:{'mass_number':46,'atomicWeight':45.99957,'abundance':0,'nuclearSpin':0},
47:{'mass_number':47,'atomicWeight':47.00762,'abundance':0,'nuclearSpin':0},
48:{'mass_number':48,'atomicWeight':48.01299,'abundance':0,'nuclearSpin':0},
49:{'mass_number':49,'atomicWeight':49.02201,'abundance':0,'nuclearSpin':0},
'most_abundant':32},

'Cl':

{'name':'Chlorine','symbol':'Cl','atomicNumber':17,'mass':35.453,
'melting_point':172.17,'boiling_point':239.18,'density':2.95,
'electron_affinity':-349,'ionization_energy_1':1252.61,
'covalent_radius':0.994,'atomic_radius':1.00,'vanderwaals_radius':1.75,'electronegativity':3.16},

'isotopes_Cl':{
28:{'mass_number':28,'atomicWeight':28.02851,'abundance':0,'nuclearSpin':0},
29:{'mass_number':29,'atomicWeight':29.01411,'abundance':0,'nuclearSpin':0},
30:{'mass_number':30,'atomicWeight':30.00477,'abundance':0,'nuclearSpin':0},
31:{'mass_number':31,'atomicWeight':30.992420,'abundance':0,'nuclearSpin':0},
32:{'mass_number':32,'atomicWeight':31.985689,'abundance':0,'nuclearSpin':0},
33:{'mass_number':33,'atomicWeight':32.9774518,'abundance':0,'nuclearSpin':0},
34:{'mass_number':34,'atomicWeight':33.97376197,'abundance':0,'nuclearSpin':0},
35:{'mass_number':35,'atomicWeight':34.96885271,'abundance':75.78,'nuclearSpin':0},
36:{'mass_number':36,'atomicWeight':35.96830695,'abundance':0,'nuclearSpin':0},
37:{'mass_number':37,'atomicWeight':36.96590260,'abundance':24.22,'nuclearSpin':0},
38:{'mass_number':38,'atomicWeight':37.96801055,'abundance':0,'nuclearSpin':0},
39:{'mass_number':39,'atomicWeight':38.9680077,'abundance':0,'nuclearSpin':0},
40:{'mass_number':40,'atomicWeight':39.970420,'abundance':0,'nuclearSpin':0},
41:{'mass_number':41,'atomicWeight':40.970650,'abundance':0,'nuclearSpin':0},
42:{'mass_number':42,'atomicWeight':41.97317,'abundance':0,'nuclearSpin':0},
43:{'mass_number':43,'atomicWeight':42.97420,'abundance':0,'nuclearSpin':0},
44:{'mass_number':44,'atomicWeight':43.97854,'abundance':0,'nuclearSpin':0},
45:{'mass_number':45,'atomicWeight':44.97970,'abundance':0,'nuclearSpin':0},
46:{'mass_number':46,'atomicWeight':45.98412,'abundance':0,'nuclearSpin':0},
47:{'mass_number':47,'atomicWeight':46.98795,'abundance':0,'nuclearSpin':0},
48:{'mass_number':48,'atomicWeight':47.99485,'abundance':0,'nuclearSpin':0},
49:{'mass_number':49,'atomicWeight':48.99989,'abundance':0,'nuclearSpin':0},
50:{'mass_number':50,'atomicWeight':50.00773,'abundance':0,'nuclearSpin':0},
51:{'mass_number':51,'atomicWeight':51.01353,'abundance':0,'nuclearSpin':0},
'most_abundant':35},

'Ar':

{'name':'Argon','symbol':'Ar','atomicNumber':18,'mass':39.948,
'melting_point':83.95,'boiling_point':87.45,'density':1.66,
'electron_affinity':29.0,'ionization_energy_1':1522.3194,
'covalent_radius':1.91,'atomic_radius':0.71,'vanderwaals_radius':1.88,'electronegativity':0},

'isotopes_Ar':{
36:{'mass_number':36,'atomicWeight':35.96754628,'abundance':0.3365,'nuclearSpin':0},
38:{'mass_number':38,'atomicWeight':37.9627322,'abundance':0.0632,'nuclearSpin':0},
40:{'mass_number':40,'atomicWeight':39.962383123,'abundance':99.6003,'nuclearSpin':0},
'most_abundant':39},

'K':

{'name':'Potassium','symbol':'K','atomicNumber':19,'mass':39.0983,
'melting_point':336.8,'boiling_point':1033,'density':0.86,
'electron_affinity':0,'ionization_energy_1':0,
'covalent_radius':2.03,'atomic_radius':2.20,'vanderwaals_radius':2.75,'electronegativity':0.82},

'isotopes_K':{
39:{'mass_number':39,'atomicWeight':38.9637,'abundance':93.2581,'nuclearSpin':0},
'most_abundant':39},

'Ca':

{'name':'Calcium','symbol':'Ca','atomicNumber':20,'mass':40.078,
'melting_point':1112,'boiling_point':1757,'density':1.74,
'electron_affinity':10.0,'ionization_energy_1':590.5158,
'covalent_radius':1.74,'atomic_radius':1.80,'vanderwaals_radius':0,'electronegativity':1.0},

'isotopes_Ca':{
40:{'mass_number':40,'atomicWeight':39.9625912,'abundance':96.941,'nuclearSpin':0},
42:{'mass_number':42,'atomicWeight':41.9586183,'abundance':0.647,'nuclearSpin':0},
43:{'mass_number':43,'atomicWeight':42.9587668,'abundance':0.135,'nuclearSpin':0},
44:{'mass_number':44,'atomicWeight':43.9554811,'abundance':2.086,'nuclearSpin':0},
46:{'mass_number':46,'atomicWeight':45.9536928,'abundance':0.004,'nuclearSpin':0},
48:{'mass_number':48,'atomicWeight':46.9545465,'abundance':0.187,'nuclearSpin':0},
'most_abundant':40},

'Sc':

{'name':'Scandium','symbol':'Sc','atomicNumber':21,'mass':44.955910,
'melting_point':1814,'boiling_point':3109,'density':2.99,
'electron_affinity':-18,'ionization_energy_1':631.764,
'covalent_radius':1.44,'atomic_radius':1.606,'vanderwaals_radius':0,'electronegativity':1.36},

'isotopes_Sc':{
45:{'mass_number':45,'atomicWeight':44.9559102,'abundance':100,'nuclearSpin':0},
'most_abundant':45},

'Br':

{'name':'Bromine','symbol':'Br','atomicNumber':35,'mass':79.904,
'melting_point':265.95,'boiling_point':331.85,'density':3.14,
'electron_affinity':-325,'ionization_energy_1':1141.23,
'covalent_radius':1.145,'atomic_radius':1.15,'vanderwaals_radius':1.85,'electronegativity':2.96},

'isotopes_Br':{
67:{'mass_number':67,'atomicWeight':66.964795,'abundance':0,'nuclearSpin':0},
68:{'mass_number':68,'atomicWeight':67.9582558,'abundance':0,'nuclearSpin':0},
69:{'mass_number':69,'atomicWeight':68.9501834,'abundance':0,'nuclearSpin':0},
70:{'mass_number':70,'atomicWeight':69.9446239,'abundance':0,'nuclearSpin':0},
71:{'mass_number':71,'atomicWeight':70.9392532,'abundance':0,'nuclearSpin':0},
72:{'mass_number':72,'atomicWeight':71.9365028,'abundance':0,'nuclearSpin':0},
73:{'mass_number':73,'atomicWeight':72.9317914,'abundance':0,'nuclearSpin':0},
74:{'mass_number':74,'atomicWeight':73.92989116,'abundance':0,'nuclearSpin':0},
75:{'mass_number':75,'atomicWeight':74.92577615,'abundance':0,'nuclearSpin':0},
76:{'mass_number':76,'atomicWeight':75.92454210,'abundance':0,'nuclearSpin':0},
77:{'mass_number':77,'atomicWeight':76.9213803,'abundance':0,'nuclearSpin':0},
78:{'mass_number':78,'atomicWeight':77.9211464,'abundance':0,'nuclearSpin':0},
79:{'mass_number':79,'atomicWeight':78.918337620,'abundance':50.697,'nuclearSpin':0},
80:{'mass_number':80,'atomicWeight':79.918530020,'abundance':0,'nuclearSpin':0},
81:{'mass_number':81,'atomicWeight':80.9162913,'abundance':49.317,'nuclearSpin':0},
82:{'mass_number':82,'atomicWeight':81.9168053,'abundance':0,'nuclearSpin':0},
83:{'mass_number':83,'atomicWeight':82.9151805,'abundance':0,'nuclearSpin':0},
84:{'mass_number':84,'atomicWeight':83.91650427,'abundance':0,'nuclearSpin':0},
85:{'mass_number':85,'atomicWeight':84.91560821,'abundance':0,'nuclearSpin':0},
86:{'mass_number':86,'atomicWeight':85.91879712,'abundance':0,'nuclearSpin':0},
87:{'mass_number':87,'atomicWeight':86.92071119,'abundance':0,'nuclearSpin':0},
88:{'mass_number':88,'atomicWeight':87.92407040,'abundance':0,'nuclearSpin':0},
89:{'mass_number':89,'atomicWeight':88.92639060,'abundance':0,'nuclearSpin':0},
90:{'mass_number':90,'atomicWeight':89.93063080,'abundance':0,'nuclearSpin':0},
91:{'mass_number':91,'atomicWeight':90.93397080,'abundance':0,'nuclearSpin':0},
92:{'mass_number':92,'atomicWeight':91.93926050,'abundance':0,'nuclearSpin':0},
93:{'mass_number':93,'atomicWeight':92.9431032,'abundance':0,'nuclearSpin':0},
94:{'mass_number':94,'atomicWeight':93.9486843,'abundance':0,'nuclearSpin':0},
'most_abundant':79}
}

ElementaryParticlesDatabase = {
'e-': {'name':'electron','symbol':'e-','category':'lepton','charge':-1.0,'mass':1.0,'spin':0.5},
'e+': {'name':'positron','symbol':'e+','category':'lepton','charge':1.0,'mass':1.0,'spin':0.5},
'u-': {'name':'muon','symbol':'u-','category':'lepton','charge':-1.0,'mass':206.76828235,'spin':0.5},
't-': {'name':'tau','symbol':'t-','category':'lepton','charge':-1.0,'mass':3477.48,'spin':0.5},
'p': {'name':'proton','symbol':'p','category':'baryon','charge':1.0,'mass':1836.15267247,'spin':0.5},
'n': {'name':'neutron','symbol':'n','category':'baryon','charge':0.0,'mass':1838.6836605,'spin':0.5},
'user': {'name':'user','symbol':' ','category':' ','charge':0.,'mass':0.0,'spin':1.0}
}

UnitsDatabase = {
'Angstroms':0.52917721092 , 'Bohr': 1.889725989
}
