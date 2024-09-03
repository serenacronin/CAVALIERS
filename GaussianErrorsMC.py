import sys
sys.path.append('../astro_tools')
from gauss_tools import one_gaussian
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def mc(mapp):

    mean_vels = []
    i = 0
    while i < 1000:
        vels = []
        j = 0
        while j < 100:
            # generate a random x and y
            x = np.random.randint(low=0, high=mapp.shape[1])  # high is exclusive
            y = np.random.randint(low=0, high=mapp.shape[0])

            if np.isfinite(mapp[y,x]):
                vels.append(mapp[y,x])  # save the velocity
                j+=1
            else:  # skip if nan
                continue

        mean_vels.append(np.mean(vels))  # take the mean of the velocities
        i+=1

    return(mean_vels)

vel_outb_map = fits.open('../ngc253/muse/data/output_maps/ngc253_muse_se_niib_vel_outflow_blue.fits')[0].data
vel_outr_map = fits.open('../ngc253/muse/data/output_maps/ngc253_muse_se_niib_vel_outflow_red.fits')[0].data
vel_disk_map = fits.open('../ngc253/muse/data/output_maps/ngc253_muse_se_niib_vel_disk.fits')[0].data

mean_vels_outb = mc(vel_outb_map)
# mean_vels_outr = mc(vel_outr_map)
# mean_vels_disk = mc(vel_disk_map)

n_outb, bins_outb, _ = plt.hist(mean_vels_outb, bins=100, color= '#66CCEE', alpha=0.6, label='Outflow B')
popt_outb, pcov_outb = curve_fit(one_gaussian, bins_outb[-1], n_outb, [30, -160, 25])
print(bins_outb)
outb_range = np.linspace(bins_outb[0],bins_outb[-1],100)
plt.plot(one_gaussian(outb_range, *popt_outb),color='#66CCEE')
         
# n_disk, bins_disk, _ = plt.hist(mean_vels_disk, bins=100, color= '#AA3377', alpha=0.6, label='Disk')
# n_outr, bins_outr, _ = plt.hist(mean_vels_outr, bins=100, color= '#CCBB44', alpha=0.6, label='Outflow R')
plt.xlabel('Mean of 100 Random Velocities')
plt.ylabel('Count')
plt.show()
