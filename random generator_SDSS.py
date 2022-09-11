#This code computes random catalogue for galaxies in the SDSS data!
#For creating ra,dec it used randomSDSS package! For creating redshift it uses n(z) of original data!


#ra, dec
import randomsdss
from scipy.stats import gaussian_kde
ra, dec = randomsdss.sky_random(dr="DR16", catalog="SDSS", size=500000)

def randz(z, randcatsize):
   

    pdf = gaussian_kde(z)
    z_rand = pdf.resample(size= randcatsize)
    z_rand = randzv.reshape(-1,1)
    z_rand = randzv.reshape(-1)

    return z_rand

