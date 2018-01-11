from nbodykit.lab import *
from nbodykit.cosmology import WMAP9, LinearPower, Cosmology
import numpy

MYPlanck = Cosmology(m_ncdm=[],
               Omega0_b=0.0223/0.6774**2,
               Omega0_cdm=0.1188/0.6774**2, 
               h=0.6774)\
            .match(sigma8=0.8159)


pklin0 = LinearPower(MYPlanck, redshift=0.0)

k = numpy.logspace(-3, 2, 10000, endpoint=True)

numpy.savetxt('myplanck-z0.txt', list(zip(k, pklin0(k))))

