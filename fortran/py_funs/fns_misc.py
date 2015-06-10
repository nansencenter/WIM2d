def SDF_Bretschneider(omega,T_peak,H_sig):
   import numpy as np

   PI    = np.pi
   T     = 2*PI/omega
   om_m  = 2*PI/T_peak
   ##
   f1 = 5.0/16.0*pow(H_sig,2)*pow(om_m,4)
   f2 = 1.0/pow(omega,5)
   f3 = exp(-1.25*pow(T/T_peak,4))
   ##
   S  = f1*f2*f3

   return S
