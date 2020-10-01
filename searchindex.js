Search.setIndex({docnames:["API","data/constants","data/data","data/databases","functional/functional","functional/interspecies_correlation","functional/xc_interface","grids/angular","grids/atomic","grids/becke","grids/cubic_spline","grids/extrapolation","grids/grids","grids/lebedev","grids/multi_grid","grids/poisson_solver","grids/radial","grids/radial_transform","gto/basis_set","gto/contracted_gaussian","gto/gto","gto/primitive_gaussian","index","scf/convergence","scf/scf","scf/scf_module","solver/dft_solver","solver/hf_solver","solver/solver","system/atomic_element","system/elementary_particle","system/input_parser","system/molecular_system","system/napmo_system","system/system","system/timer","utilities/cell","utilities/int1d","utilities/ode2","utilities/utilities","wavefunction/ndpsi","wavefunction/nkinetic","wavefunction/nnuclear","wavefunction/ntwobody","wavefunction/psi_analytic","wavefunction/psi_base","wavefunction/psi_hybrid","wavefunction/psi_numeric","wavefunction/psi_optimization","wavefunction/wavefunction"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":3,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":2,"sphinx.domains.rst":2,"sphinx.domains.std":1,sphinx:56},filenames:["API.rst","data/constants.rst","data/data.rst","data/databases.rst","functional/functional.rst","functional/interspecies_correlation.rst","functional/xc_interface.rst","grids/angular.rst","grids/atomic.rst","grids/becke.rst","grids/cubic_spline.rst","grids/extrapolation.rst","grids/grids.rst","grids/lebedev.rst","grids/multi_grid.rst","grids/poisson_solver.rst","grids/radial.rst","grids/radial_transform.rst","gto/basis_set.rst","gto/contracted_gaussian.rst","gto/gto.rst","gto/primitive_gaussian.rst","index.rst","scf/convergence.rst","scf/scf.rst","scf/scf_module.rst","solver/dft_solver.rst","solver/hf_solver.rst","solver/solver.rst","system/atomic_element.rst","system/elementary_particle.rst","system/input_parser.rst","system/molecular_system.rst","system/napmo_system.rst","system/system.rst","system/timer.rst","utilities/cell.rst","utilities/int1d.rst","utilities/ode2.rst","utilities/utilities.rst","wavefunction/ndpsi.rst","wavefunction/nkinetic.rst","wavefunction/nnuclear.rst","wavefunction/ntwobody.rst","wavefunction/psi_analytic.rst","wavefunction/psi_base.rst","wavefunction/psi_hybrid.rst","wavefunction/psi_numeric.rst","wavefunction/psi_optimization.rst","wavefunction/wavefunction.rst"],objects:{"":{angular:[7,0,0,"-"],atomic:[8,0,0,"-"],atomic_element:[29,0,0,"-"],basis_set:[18,0,0,"-"],becke:[9,0,0,"-"],cell:[36,0,0,"-"],constants:[1,0,0,"-"],contracted_gaussian:[19,0,0,"-"],convergence:[23,0,0,"-"],cubic_spline:[10,0,0,"-"],databases:[3,0,0,"-"],dft_solver:[26,0,0,"-"],elementary_particle:[30,0,0,"-"],extrapolation:[11,0,0,"-"],hf_solver:[27,0,0,"-"],input_parser:[31,0,0,"-"],int1d:[37,0,0,"-"],interspecies_correlation:[5,0,0,"-"],lebedev:[13,0,0,"-"],molecular_system:[32,0,0,"-"],multi_grid:[14,0,0,"-"],napmo_system:[33,0,0,"-"],ndpsi:[40,0,0,"-"],nkinetic:[41,0,0,"-"],nnuclear:[42,0,0,"-"],ntwobody:[43,0,0,"-"],ode2:[38,0,0,"-"],poisson_solver:[15,0,0,"-"],primitive_gaussian:[21,0,0,"-"],psi_analytic:[44,0,0,"-"],psi_base:[45,0,0,"-"],psi_hybrid:[46,0,0,"-"],psi_numeric:[47,0,0,"-"],psi_optimization:[48,0,0,"-"],radial:[16,0,0,"-"],radial_transform:[17,0,0,"-"],scf:[25,0,0,"-"],timer:[35,0,0,"-"],xc_interface:[6,0,0,"-"]},"angular.AngularGrid":{integrate:[7,2,1,""],lorder:[7,2,1,""],spherical:[7,2,1,""]},"atomic.AtomicGrid":{evaluate_expansion:[8,2,1,""],integrate:[8,2,1,""],spherical_average:[8,2,1,""],spherical_expansion:[8,2,1,""],symbol:[8,2,1,""]},"atomic_element.AtomicElement":{is_quantum:[29,2,1,""]},"basis_set.BasisSet":{compute:[18,2,1,""],deriv:[18,2,1,""],load_file:[18,2,1,""],update:[18,2,1,""]},"becke.BeckeGrid":{becke_weights:[9,2,1,""],evaluate_decomposition:[9,2,1,""],integrate:[9,2,1,""],ncenter:[9,2,1,""],show:[9,2,1,""]},"contracted_gaussian.ContractedGaussian":{compute:[19,2,1,""],l:[19,2,1,""],origin:[19,2,1,""],overlap:[19,2,1,""]},"convergence.Convergence":{damping:[23,2,1,""]},"cubic_spline.CubicSpline":{deriv:[10,2,1,""]},"dft_solver.DFT":{compute:[26,2,1,""],compute_analytic:[26,2,1,""],compute_hybrid:[26,2,1,""],compute_numeric:[26,2,1,""],get:[26,2,1,""],pce:[26,2,1,""]},"elementary_particle.ElementaryParticle":{is_quantum:[30,2,1,""]},"extrapolation.Extrapolation":{to_string:[11,2,1,""]},"extrapolation.PotentialExtrapolation":{to_string:[11,2,1,""]},"extrapolation.PowerExtrapolation":{to_string:[11,2,1,""]},"hf_solver.HF":{compute:[27,2,1,""],compute_analytic:[27,2,1,""],compute_hybrid:[27,2,1,""],compute_numeric:[27,2,1,""],get:[27,2,1,""],pce:[27,2,1,""]},"input_parser.InputParser":{load_basis:[31,2,1,""],load_code:[31,2,1,""],load_functional:[31,2,1,""],load_molecule:[31,2,1,""],load_scf:[31,2,1,""]},"int1d.CubicIntegrator1D":{get_weights:[37,2,1,""],npoint_min:[37,4,1,""]},"int1d.Integrator1D":{get_weights:[37,2,1,""],npoint_min:[37,4,1,""]},"int1d.SimpsonIntegrator1D":{get_weights:[37,2,1,""],npoint_min:[37,4,1,""]},"int1d.StubIntegrator1D":{get_weights:[37,2,1,""],npoint_min:[37,4,1,""]},"int1d.TrapezoidIntegrator1D":{get_weights:[37,2,1,""],npoint_min:[37,4,1,""]},"molecular_system.MolecularSystem":{add_atom:[32,2,1,""],add_elementary_particle:[32,2,1,""],add_nucleus:[32,2,1,""],get_basis:[32,2,1,""],get_species:[32,2,1,""],point_charges:[32,2,1,""],set_charges:[32,2,1,""],size_particles:[32,2,1,""],size_species:[32,2,1,""]},"multi_grid.MultiGrid":{add_grid:[14,2,1,""],get_common_points:[14,2,1,""],get_grid:[14,2,1,""],get_grid_id:[14,2,1,""],ngrids:[14,2,1,""],nspecies:[14,2,1,""],show:[14,2,1,""]},"napmo_system.NAPMO":{build_system:[33,2,1,""],exec_code:[33,2,1,""],solve:[33,2,1,""]},"primitive_gaussian.PrimitiveGaussian":{compute:[21,2,1,""],l:[21,2,1,""],origin:[21,2,1,""],overlap:[21,2,1,""]},"psi_analytic.PSIA":{build_fock:[44,2,1,""],compute_2body:[44,2,1,""],compute_c_2species_grid:[44,2,1,""],compute_coupling:[44,2,1,""],compute_guess:[44,2,1,""],compute_hcore:[44,2,1,""],compute_kinetic:[44,2,1,""],compute_nuclear:[44,2,1,""],compute_overlap:[44,2,1,""],compute_transformation:[44,2,1,""],compute_xc_grid:[44,2,1,""],compute_xc_matrix:[44,2,1,""]},"psi_base.PSIB":{nbasis:[45,2,1,""],pce:[45,2,1,""],plot_obj:[45,2,1,""],sid:[45,2,1,""],symbol:[45,2,1,""]},"psi_hybrid.PSIH":{compute_2body:[46,2,1,""],compute_coupling:[46,2,1,""]},"psi_numeric.PSIN":{build_fock:[47,2,1,""],compute_1body:[47,2,1,""],compute_2body:[47,2,1,""],compute_coupling:[47,2,1,""],compute_guess:[47,2,1,""],compute_hcore:[47,2,1,""],compute_kinetic:[47,2,1,""],compute_nuclear:[47,2,1,""],compute_overlap:[47,2,1,""],normalize:[47,2,1,""],optimize_psi:[47,2,1,""]},"psi_optimization.PSIO":{build_fock:[48,2,1,""],compute_1body:[48,2,1,""],compute_2body:[48,2,1,""],compute_cor2species:[48,2,1,""],compute_coupling:[48,2,1,""],compute_guess:[48,2,1,""],compute_hcore:[48,2,1,""],compute_kinetic:[48,2,1,""],compute_nuclear:[48,2,1,""],compute_overlap:[48,2,1,""],normalize:[48,2,1,""],optimize:[48,2,1,""]},"radial.RadialGrid":{integrate:[16,2,1,""]},"scf.SCF":{compute_energy:[25,2,1,""],compute_energy_components:[25,2,1,""],compute_energy_single:[25,2,1,""],energy:[25,2,1,""],get:[25,2,1,""],hmulti:[25,2,1,""],multi:[25,2,1,""],nmulti:[25,2,1,""],nsingle:[25,2,1,""],pce:[25,2,1,""],show_results:[25,2,1,""],single:[25,2,1,""]},"timer.Timer":{show_block:[35,2,1,""],show_summary:[35,2,1,""],timeblock:[35,2,1,""]},"xc_interface.Functional":{compute_correlation:[6,2,1,""],compute_exchange:[6,2,1,""],show:[6,2,1,""]},angular:{AngularGrid:[7,1,1,""]},atomic:{AtomicGrid:[8,1,1,""]},atomic_element:{AtomicElement:[29,1,1,""]},basis_set:{BasisSet:[18,1,1,""]},becke:{BeckeGrid:[9,1,1,""]},cell:{Cell:[36,1,1,""]},contracted_gaussian:{ContractedGaussian:[19,1,1,""]},convergence:{Convergence:[23,1,1,""]},cubic_spline:{CubicSpline:[10,1,1,""]},databases:{AtomicElementsDatabase:[3,3,1,""],CouplingConstantsDatabase:[3,3,1,""],ElementaryParticlesDatabase:[3,3,1,""]},dft_solver:{DFT:[26,1,1,""]},elementary_particle:{ElementaryParticle:[30,1,1,""]},extrapolation:{CuspExtrapolation:[11,1,1,""],Extrapolation:[11,1,1,""],PotentialExtrapolation:[11,1,1,""],PowerExtrapolation:[11,1,1,""]},hf_solver:{HF:[27,1,1,""]},input_parser:{InputParser:[31,1,1,""],extract_keywork:[31,3,1,""],fix_casting:[31,3,1,""]},int1d:{CubicIntegrator1D:[37,1,1,""],Integrator1D:[37,1,1,""],SimpsonIntegrator1D:[37,1,1,""],StubIntegrator1D:[37,1,1,""],TrapezoidIntegrator1D:[37,1,1,""]},interspecies_correlation:{epc17_2:[5,3,1,""],isc_functional_selector:[5,3,1,""]},lebedev:{lebedev_get_order:[13,3,1,""]},molecular_system:{MolecularSystem:[32,1,1,""]},multi_grid:{MultiGrid:[14,1,1,""]},napmo_system:{NAPMO:[33,1,1,""]},ndpsi:{compute_dpsi:[40,3,1,""]},nkinetic:{compute_kinetic:[41,3,1,""]},nnuclear:{compute_nuclear:[42,3,1,""]},ntwobody:{compute_coulomb:[43,3,1,""]},ode2:{solve_ode2:[38,3,1,""]},poisson_solver:{poisson_solver:[15,3,1,""]},primitive_gaussian:{PrimitiveGaussian:[21,1,1,""]},psi_analytic:{PSIA:[44,1,1,""]},psi_base:{PSIB:[45,1,1,""]},psi_hybrid:{PSIH:[46,1,1,""]},psi_numeric:{PSIN:[47,1,1,""]},psi_optimization:{PSIO:[48,1,1,""]},radial:{RadialGrid:[16,1,1,""]},radial_transform:{ChebyshevRadialTransform:[17,1,1,""],IndentityRadialTransform:[17,1,1,""],PowerRadialTransform:[17,1,1,""],RadialTransform:[17,1,1,""]},scf:{SCF:[25,1,1,""]},timer:{Timer:[35,1,1,""]},xc_interface:{Functional:[6,1,1,""]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","function","Python function"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:function","4":"py:attribute"},terms:{"00794":3,"084":3,"1021":5,"1039":3,"110":9,"1312":3,"15267247":1,"15432897":18,"1836":1,"1838":1,"1986":21,"1988":[9,15,17],"2008":3,"2015":22,"2018":6,"2547":[9,17],"2832":3,"42525091":18,"52917721092":1,"5858151015155881":3,"6836605":1,"7b01442":5,"889725989":1,"abstract":[17,30],"case":[8,14,18,25],"class":[6,7,8,9,10,11,14,16,17,18,19,21,23,25,26,27,29,30,31,32,33,35,36,37,44,45,46,47,48],"default":[6,9,10,22,25,26,27,29,32],"float":[8,17,40],"function":[0,5,6,7,8,9,10,11,17,18,19,21,22,28,31,38,40,41,43],"import":[3,14,31],"int":[8,9,11,13,14,15,17,31,32,40,41,43],"new":32,"return":[5,7,8,9,10,11,13,14,15,18,25,26,27,32,37,38,41,42,43],"true":[3,25,26,27,31,33],"while":14,And:1,For:17,Not:17,ODE:38,One:14,RKS:22,The:[9,10,11,14,17,19,21,22,25,26,27,36,38,42,43,45],These:[17,38],UKS:22,_den:15,a_i:[19,21],a_x:[19,21],a_z:[19,21],about:[14,35],accordingli:32,account:43,acs:5,action:41,actual:17,add:[14,32],add_atom:32,add_elementary_particl:32,add_grid:14,add_nucleu:32,added:[14,32],addit:31,adjac:11,algorithm:[23,37],all:[3,8,9,18,25,29,35,42,47,48],along:8,alpha:[17,32],alreadi:19,also:19,analyt:[21,44,47,48],angstrom:[1,29,30,32],angstrom_to_bohr:1,angular:[0,9,11,12,18,19,21],angulargrid:7,ani:43,api:22,apmo:[22,26,27],appreci:22,approach:22,approxim:[8,9],arg:[7,8,48],argument:[10,14],arrai:[7,8,9,10,14,15,18,19,21,29,30],atom:[0,3,9,12,17,29,30,32],atomic_el:[0,34],atomic_numb:[3,42],atomic_radii:3,atomic_radii_2:3,atomic_symbol:[8,16],atomicel:29,atomicelementsdatabas:3,atomicgrid:8,attract:42,author:22,avail:17,averag:8,axi:17,b801115j:3,base:[3,7,17,22,37],basi:[18,22,31,32,45,47],basis_data:18,basis_fil:[18,32],basis_nam:[18,32],basis_set:[0,20,32],basisset:[18,32],bcs:38,beck:[0,12,15,17,40,43,47,48],becke_weight:9,beckegrid:[9,14,15,40,41,42,43],beg:31,being:[29,30],besid:31,beta:32,beta_psi:44,between:[1,5,14,19,21],block:33,bodi:[44,46,47,48],bohr:[1,29,30,32],bohr_to_angstrom:1,boiling_point:3,bool:[25,26,27,31,32,44,46,47,48],both:14,boundari:[11,36,38],bse:18,build:[33,43,44,47,48],build_fock:[44,47,48],build_neumann:38,build_system:33,c_i:19,calcul:[5,7,8,9,11,14,15,19,21,25,26,27,31,33,35,41,43,44,45,46,47,48],call:17,can:[11,22,47],cartesian:[14,19,21,29,31,32],categori:3,cell:[0,9,39],cell_funct:9,center:[8,17,19,21,43],certain:35,charg:[3,25,26,27,32,40,42,45],chebyshevradialtransform:17,check:22,chem:[9,15,17],choos:5,code:[22,31,33],coeffici:[19,21],combin:19,common:[5,14],compil:22,compon:[16,19,21],composit:37,comprehens:6,comput:[6,8,9,11,18,19,21,25,26,27,40,41,42,44,45,46,47,48],compute_1bodi:[47,48],compute_2bodi:[44,46,47,48],compute_analyt:[26,27],compute_c_2species_grid:44,compute_cor2speci:48,compute_correl:6,compute_coulomb:43,compute_coupl:[44,46,47,48],compute_dpsi:40,compute_energi:25,compute_energy_compon:25,compute_energy_singl:25,compute_exchang:6,compute_guess:[44,47,48],compute_hcor:[44,47,48],compute_hybrid:[26,27],compute_kinet:[41,44,47,48],compute_nuclear:[42,44,47,48],compute_numer:[26,27],compute_overlap:[44,47,48],compute_transform:44,compute_xc_grid:44,compute_xc_matrix:44,condit:[11,36,38],conrad:6,consist:[0,11,22],constant:[0,2,3,19,29,30],cont:18,contain:[1,3,22,32],context:35,continu:[10,11],contract:19,contracted_gaussian:[0,20],contractedgaussian:19,convent:[26,27],converg:[0,24],convert:31,coord:[18,19,21],coordin:[7,14,16,18,19,21,31,32],cordero:3,correa:22,correl:[6,44,48],correspond:[14,31],could:22,coulomb:43,coupl:[3,44,46,47,48],couplingconstantsdatabas:3,coval:[3,17],cover:22,coverag:22,covert:1,creat:[9,14,30],cubic:[10,37,38,43],cubic_splin:[0,9,12],cubicintegrator1d:37,cubicsplin:[9,10,38],cuda:22,cuspextrapol:[10,11],dalton:3,damp:23,data:[0,3,18,22,31,32,33],databas:[0,1,2,29,30],debug:47,decompos:43,decomposit:9,defin:[7,8,17,19,21,32,44,45,46],definit:17,delta:[40,47,48],demon2k:18,den:[15,43],denot:21,densiti:[0,3,5,6,15,28,43,44,45,47,48],densitii:6,depend:[26,27],deriv:[10,18],describ:9,desir:14,detail:[25,38,43],determin:10,deuterium:32,develop:6,dfrac:[17,19,41,42],dft:[0,28,31],dft_solver:26,dickson:15,dict:[6,9,18,25,30,32,42],differ:[17,19,32],dii:25,dim:48,dimension:36,dimes:5,direct:[44,46,47,48],directli:17,distribut:17,doctest:22,doe:18,doi:[3,5,40],don:36,doubl:[7,21,25],dtype:[19,21],each:[7,8,15,32,45],edu:22,edwin:22,efposadac:22,either:[18,47],electron:[3,5,6,21,30,32,44,47,48],electron_affin:3,electron_charg:1,electron_mass:1,electroneg:3,electrostat:40,element:[3,15,29,30,38],elementari:[30,32],elementary_particl:[0,34],elementaryparticl:30,elementaryparticlesdatabas:3,ell:[8,9,15],ell_:[8,9,15],energi:[6,25,26,27,40,45],entir:9,epc17:[5,22],epc17_2:5,equat:[11,15,43],eri:[44,46,47,48],eta:3,etc:[3,30,32],evalu:[8,9,10,18,40],evaluate_decomposit:9,evaluate_expans:8,exampl:[18,31],exchang:[6,44,48],exec_cod:33,execut:[33,35],exist:32,exp:[11,19,21],expand:[15,47],expans:[8,9,15,40,41,43],expon:[19,21],extens:6,extract:31,extract_keywork:31,extrapol:[0,5,10,12,38],factor:1,fals:[3,7,25,32,36,44,46,47,48],fernando:22,field:[0,22],file:[9,18,31,33],finit:[15,38],first:18,fix_cast:31,flavor:22,fly:[44,46,47,48],fock:[0,22,28,44,45,46,47,48],follow:[15,17,18,21,43],form:[11,19],format:18,found:22,free:22,from:[3,13,17,18,19,29,30,31,33,40,45,47],further:18,gaussian:[0,19,21,22],get:[25,26,27],get_basi:32,get_common_point:14,get_grid:14,get_grid_id:14,get_speci:32,get_weight:37,gid:14,github:22,given:[5,8,9,10,11,14,17,18,32,38],global:35,gov:[3,18,29,30],grid:[0,5,6,7,8,9,10,11,14,15,16,17,22,37,40,41,42,43,44,45,46,47,48],group:31,gsl:22,gto:[0,19,21,22,26,27],guess:[44,47,48],h_2:32,hamm:5,handl:[14,25,29],harmon:43,hartre:[0,22,28,44,45,46],has:[14,18],have:[5,10,14,19,38,42],hcore:[44,47,48],help:22,henc:19,hf_solver:27,hmulti:25,http:[3,5,18,22,29,30],hybrid:[26,27,46],hydrogen:3,ident:[10,38],implement:[17,22,23,25],includ:[19,45],indentityradialtransform:17,index:[9,14,22],indic:19,inform:[3,6,9,14,18,29,30,31,32,35,42],initi:47,initvoid:36,input:[31,33,38],input_pars:[0,34],inputpars:[31,33],instal:22,int1d:[0,39],int32:[19,21],int_:8,integ:21,integr:[7,8,9,16,17,19,21,37,44,46,47,48],integrator1d:37,interfac:[11,18,30],interpol:43,interspeci:5,interspecies_correl:[0,4],interv:10,invers:40,ionization_energy_1:3,is_quantum:[3,29,30],isc_functional_selector:5,item:31,iter:[18,25],its:31,jpclett:5,just:[19,31],kappa:3,kei:[18,25,26,27,42],keyword:31,kind:32,kinet:[41,44,47,48],know:14,kwarg:8,label:[35,45],lack:18,lambda:3,last:11,lda:22,lebedev:[0,7,12],lebedev_get_ord:13,left:11,lehtola:6,length:45,lepton:[3,30],libint:22,librari:[6,22],libxc:[6,22],limit:40,linear:[19,37,48],list:[9,18,32,42],lmax:[8,15,40,41,43],load:[18,31],load_basi:31,load_cod:31,load_fil:18,load_funct:31,load_molecul:31,load_scf:31,local:22,logic:8,lorder:7,lsize:[8,15],make:22,manag:[33,35],marker:45,marqu:6,mass:[3,40],mass_inv:40,matric:47,matrix:[44,45,46,47,48],max:[8,9,15],maximum:[8,15,43],mean:17,melting_point:3,method:[18,22,44],mhpc:22,micael:6,miguel:6,minim:47,modul:22,molecul:[8,9,15,17,31,32],molecular:[9,14,15,17,32,40,41,42,43],molecular_system:[0,34],molecularsystem:[9,32,33],moment:[19,21],momentum:[11,19,21],more:22,multi:[14,17,25,43],multi_grid:[0,12],multicent:9,multigrid:14,multipl:[14,32],muon:[3,30,32],must:[5,10,14,38],n_angular:9,n_radial:9,n_x:[19,21],n_y:[19,21],n_z:[19,21],nabla:41,name:[3,5],nang:8,napmo:[0,3,31,32,34,43],napmo_system:33,natur:43,nbasi:45,ncenter:9,ndarai:41,ndarrai:[5,6,7,8,9,10,14,15,18,19,21,29,30,32,40,41,42,43,47],ndim:[45,47],ndpsi:[0,49],need:36,neg:21,neo:22,neutron_charg:1,neutron_mass:1,new_dx:10,new_x:10,ngrid:14,nist:[3,29,30],nkinet:[0,49],nmulti:25,nnuclear:[0,49],non:21,none:[8,9,10,15,16,18,25,26,27,31,32,36,37,38,44,45,47,48],normal:[19,47,48],nosetest:22,note:[5,14,22],now:31,npoint:[17,37],npoint_min:37,nrad:[8,15],nsingl:25,nspeci:14,ntwobodi:[0,49],nuclear:[42,44,47,48],nuclei:[5,29],nucleu:32,number:[8,9,13,14,17,21,26,27,32],numer:[9,15,17,25,26,27,44,47],numpi:10,obara:21,obj:45,object:[6,9,10,11,14,19,21,25,26,27,32,33,38,44,45,46,47,48],obtain:[15,47],ode2:[0,39],oliveira:6,omega:8,omn:22,omp:22,one:[8,22,25],onli:[17,31,44],open_shel:32,oper:[41,42,44,45,46],optim:48,optimize_psi:47,optio:6,option:[6,9,10,25,26,27,29,30,31,32,38,44,45],orbit:[0,19,21,40,41,45,47,48],order:[7,8,13,15,38],ordinari:37,org:[3,5],origin:[8,18,19,21,29,30,32,42],orthogon:36,other:[1,19,21,44,46,47,48],other_psi:[25,44,46,47,48],other_wf:48,output:[8,9,10,14,38],outsid:10,over:[7,14,41],overlap:[19,21,44,47,48],owner:[8,45],packag:22,page:22,paper:[17,40,43,47,48],paramet:[5,6,7,8,9,10,11,13,14,15,17,18,19,21,25,26,27,29,30,31,32,33,38,40,41,42,43,44,46,47,48],pars:33,parser:31,particl:[9,18,30,32],particlesfract:3,path:18,pce:[25,26,27,45],perform:[7,8,9,16,19,21,25,33],period:36,phi:[19,21,40,47,48],phi_:[19,21,41],phy:[9,15,17],physic:[1,3,29,30],plot_obj:45,pnl:18,point:[5,7,8,9,10,11,13,14,15,17,19,25,26,27,32,42],point_charg:[32,42,44,46],poisson:[11,15,43],poisson_solv:[0,12,43],polyatom:[9,15,17],polynomi:11,portal:18,posada:22,posit:31,positron:32,potenti:[6,11,15,40,43,44,45,47,48],potentialextrapol:11,power:[11,17],powerextrapol:11,powerradialtransform:17,pprint:[25,26,27,33],prefactor:11,prerequisit:22,present:18,prim:18,primit:[19,21],primitive_gaussian:[0,20],primitivegaussian:21,print:[6,9,25,26,27,35],prioriti:14,procedur:[15,25,43],product:8,program:22,progress:[25,26,27],properli:19,properti:[7,8,9,14,19,21,25,26,27,29,30,32,45],proton_charg:1,proton_mass:1,provid:[6,14,22],psi:[25,40,41,48],psi_analyt:[0,49],psi_bas:[0,49],psi_grid:47,psi_hybrid:[0,49],psi_numer:[0,49],psi_optim:[0,49],psia:[44,46],psib:45,psih:46,psin:[47,48],psio:48,psix:47,publish:19,pylibxc:[6,22],python:[6,11,22,30],quadratur:[7,13],quantum:[21,26,27,29,30,32],quatum:3,r_0:17,r_a:42,r_i:[17,42],radial:[0,8,9,12,17],radial_transform:[0,12],radialgrid:16,radialtransform:17,radii:[3,17],radiu:17,recent:[6,22],refer:[5,6,9,15,17,21],relat:[3,6,29],relev:31,remain:11,represent:36,repuls:42,requir:17,residu:40,respresent:11,result:[19,25],revisit:3,rhf:22,rho:6,rhoe:5,rhon:5,right:11,rmax:17,rmin:17,robust:22,rtransform:[8,9,10,16,38],rule:37,run:22,rvec:36,saika:21,same:[5,10,14,19,38],scf:[0,22,26,27,31,47,48],scheme:[9,17],schiffer:5,scipi:22,search:22,second:[10,38],see:[9,18,38,43],segment:[8,16],self:[0,19,21,22],serial:22,set:[11,18,21,22,32],set_charg:32,sever:23,shape:[8,14,15,18],shell:19,should:[18,32,42],show:[6,9,14,25],show_block:35,show_result:25,show_summari:35,shown:14,sid:[32,45],side:11,simpson:37,simpsonintegrator1d:37,singl:25,size:[8,10,14,15,16,17,32,41,43],size_particl:32,size_speci:32,softwar:6,solut:[11,15,38],solv:[33,38,43,48],solve_ode2:38,solver:[0,22,38],some:[1,47],sometim:14,sourc:[15,33,43],space:17,speci:[3,6,9,14,25,26,27,29,30,32,44,45,46,47,48],specif:14,specifi:10,sph_exp:15,sphere:7,spheric:[7,8,9,15,16,40,41,43],spherical_averag:8,spherical_expans:8,spin:3,spin_boson:1,spin_down_electron:1,spin_electron:1,spin_fermion:1,spin_up_electron:1,spline:[10,11,37,38,43],split:32,sqrt:19,start:[17,31],steigemann:6,str:[5,6,14,18,29,30,31,32],string:[11,31],stubintegrator1d:37,suggest:22,suitabl:11,sum:21,sum_:[8,9,19,42],support:[22,36],susi:6,symbol:[3,6,8,14,18,29,30,32,45],symbol_a:14,symbol_b:14,system:[0,8,14,16,22,25,26,27,31,32,33,36],take:43,taken:3,task:33,templ:22,term:19,test:[17,31],than:22,thei:10,them:[5,19],theori:[0,6,28],theta:[8,9],thi:[1,3,6,9,10,14,17,18,22,25,29,43,47],time:[19,21,35],timeblock:35,timer:[0,34],to_str:11,total:[14,25],total_mass:[44,46],toward:10,transact:3,transform:[10,17],trapezoid:37,trapezoidintegrator1d:37,treat:[29,30],trial:41,two:[14,19,21,44,46,47,48],type:[0,7,8,9,13,14,15,18,19,21,22,31,32,38,41,42,43],uhf:22,unal:22,uniform:17,uniqu:14,unit:[7,29,30,32],unitsdatabas:[0,2],unnorm:21,updat:18,use:[8,17,22],used:[6,10,14],user:22,uses:[6,38],using:[7,8,9,26,27,43,44,47,48],util:0,utilitit:[0,22],v_tot:40,valid:[29,30,32],valu:[7,8,9,10,17,19,21,29,30,32],variabl:31,variat:48,varphi:[8,9],vector:[21,36],version:22,virtual:45,wai:[8,19],wavefunct:[0,22,25,44,46,47,48],weight:[9,37],welcom:22,when:10,where:[8,15,19,21],whether:[7,8,25,26,27,29,30,31,32,44,46,47,48],which:[10,11,14,18],write:21,xc_interfac:[0,4],you:5,zeta:[19,21],zeta_i:19},titles:["API","constants","Data","databases","Functional","interspecies_correlation","xc_interface","angular","atomic","becke","cubic_spline","extrapolation","Grids","lebedev","multi_grid","poisson_solver","radial","radial_transform","basis_set","contracted_gaussian","Gaussian Type Orbitals (GTO)","primitive_gaussian","nAPMO\u2019s documentation.","convergence","Self-Consistent Field (SCF)","scf","Density Functional Theory (DFT) Solver","Hartree-Fock (HF) Solver","Solver","atomic_element","elementary_particle","input_parser","molecular_system","napmo","System","timer","cell","int1d","ode2","Utilitites","ndpsi","nkinetic","nnuclear","ntwobody","psi_analytic","psi_base","psi_hybrid","psi_numeric","psi_optimization","WaveFunction"],titleterms:{"function":[4,26],angular:7,ani:22,api:0,atom:8,atomic_el:29,basis_set:18,beck:9,cell:36,consist:24,constant:1,content:22,contracted_gaussian:19,converg:23,cubic_splin:10,data:2,databas:3,densiti:26,dft:26,document:22,elementary_particl:30,extrapol:11,field:24,fock:27,gaussian:20,grid:12,gto:20,hartre:27,indic:22,input_pars:31,int1d:37,interspecies_correl:5,lebedev:13,molecular:22,molecular_system:32,multi_grid:14,napmo:[22,33],ndpsi:40,nkinet:41,nnuclear:42,ntwobodi:43,numer:22,ode2:38,orbit:[20,22],particl:22,poisson_solv:15,primitive_gaussian:21,psi_analyt:44,psi_bas:45,psi_hybrid:46,psi_numer:47,psi_optim:48,radial:16,radial_transform:17,scf:[24,25],self:24,solver:[26,27,28],system:34,tabl:22,theori:26,timer:35,type:20,unitsdatabas:1,utilitit:39,wavefunct:49,xc_interfac:6}})