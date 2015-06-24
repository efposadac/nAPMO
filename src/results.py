import matplotlib.pyplot as plt
import numpy as np
import time
import os

nprocs = np.array([i for i in range(1, 9)])
times = []

for nproc in nprocs:
    start_time = time.time()
    os.system('export OMP_NUM_THREADS='+str(nproc)+'; ./test')
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
