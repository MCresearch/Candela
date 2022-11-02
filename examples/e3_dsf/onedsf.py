import numpy as np
bplot=True
try:
    import matplotlib.pyplot as plt
except:
    bplot=False

import warnings
warnings.filterwarnings('ignore')

interval = 0.00025  # ps
tcut = 0.4  # ps
N_max = 25
#dw = 2*pi/(N_max * tcut) * 0.658212 #meV

isffile = 'isf.txt'
outfile = 'dsf.txt'

isf = np.loadtxt(isffile)

if isf[-1,0]/tcut < N_max:
    n_max = int(len(isf[:, 1]) * N_max * tcut/isf[-1,0]) - len(isf[:, 1])
    add_isf = np.zeros((n_max, 2))
    isf = np.concatenate((isf, add_isf), axis=0)

sigma =  tcut / 3
for ii in range(isf[:, 1].reshape(-1).size):
    if isf[ii][0] > tcut:
        phi = -1 * (isf[ii][0] - tcut) / sigma
        if phi < -100:
            isf[ii][1] = 0
        else:
            isf[ii][1] = isf[ii][1] * np.exp(phi)

dsf = np.real(np.fft.rfft(isf[:, 1].reshape(-1))) * interval / np.pi / 0.658212
freq = np.fft.rfftfreq(isf[:, 1].reshape(-1).size, interval)
freq = freq * 2 * np.pi * 0.658212

print('w max')
max_dsf = np.max(dsf)
max_freq = float(freq[np.where(dsf==np.max(dsf))])
print('%s %s'%(max_freq, max_dsf))


outdsf = np.dstack((freq, dsf))
np.savetxt(outfile, outdsf.reshape(-1,2))

if bplot:
    plt.plot(freq, dsf, label='DSF')
    plt.xlim(0, 1000)
    plt.legend()
    plt.show()
