import numpy as np
import matplotlib.pyplot as plt

tau0  = 40.e3  # 
mu    = .7     # 

tau_t = 42.e3  # max sigN
nu    = .3
alp   = (1+nu)/(1-nu)

sigN  = np.arange(-100.,60.,1.)*1.e3
taup  = tau0-mu*sigN
taum  = -tau0+mu*sigN

tauw  = alp*sigN
sigNp = tau0/(alp+mu)
sigNm = -tau0/(alp-mu)
sigC  = np.array([sigNp,sigNm])
tauC  = alp*sigC

sig11_p     = 2./(1-nu)*sigNp
epsc_coeff  = (1-nu*nu)*sig11_p
Yvec        = 1.e9*np.arange(.5,10.1,.1)

f  = plt.figure()

# MC envelope
ax1 = f.add_subplot(121)
ax1.plot(sigN/1.e3,taup/1.e3,'--r')
ax1.plot(sigN/1.e3,taum/1.e3,'--r')
ax1.plot(sigN/1.e3,0*sigN,'k')
ax1.plot(sigN/1.e3,tauw/1.e3,'b')
ax1.plot(sigC/1.e3,tauC/1.e3,'.r')
ax1.set_xlabel(r"$\sigma_N$, kPa")
ax1.set_ylabel(r"$\tau$, kPa")

# epsc vs Y
ax2 = f.add_subplot(122)
ax2.plot(Yvec/1.e9,(epsc_coeff/1.e-5)/Yvec,'b')
ax2.set_xlabel(r"$Y$, GPa")
ax2.set_ylabel(r"$10^5\epsilon_c$, ")

# save and close fig
f.savefig('out/MC.eps')
ax1.cla()
ax2.cla()
f.clf()
