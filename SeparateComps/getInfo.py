from astropy.io import fits # type: ignore
from astropy.wcs import wcs # type: ignore
import pandas as pd # type: ignore
import sys
sys.path.append('/Users/serenac/Desktop/research/astro_tools')
import axes_offset # type: ignore

global pixscale; pixscale = 0.2  # arcsec
global chan_width; chan_width = 1.25 # Angstrom
global num_obs; num_obs = 150
global free_params1; free_params1 = 6
global free_params2; free_params2 = 12
global free_params3; free_params3 = 18
global savepath; savepath = '/Users/serenac/Desktop/research/ngc253/muse/2024July2/separate_components/'
savepath1 = '/Users/serenac/Desktop/research/ngc253/muse/2024July2/fits1_total/'
savepath2 = '/Users/serenac/Desktop/research/ngc253/muse/2024July2/fits2_total/'
savepath3 = '/Users/serenac/Desktop/research/ngc253/muse/2024July2/fits3_total/'
infile1 = '%sfits1_reordered.txt' % savepath1
infile2 = '%sfits2_reordered.txt' % savepath2
infile3 = '%sfits3_reordered.txt' % savepath3

# read in the results of the fitting
global fits1; fits1 = pd.read_csv(infile1)
global fits2; fits2 = pd.read_csv(infile2)
global fits3; fits3 = pd.read_csv(infile3)

# get info of original data
og = '/Users/serenac/Desktop/research/ngc253/muse/data/NGC253_MUSE_SE_Fitted.fits'
hdu = fits.open(og)[1]
og_data = hdu.data
global og_hdr; og_hdr =  
global y, x; _, y, x = og_data.shape
global w; w = wcs.WCS(hdu.header, naxis=2).celestial
global offset_ra, offset_dec; offset_ra, offset_dec = axes_offset.get_offset(og_data, w)
global extent; extent = axes_offset.get_extent(og_data, x_dim=2, y_dim=1, pixscale=pixscale)

# open the BIC_PHYS test fits file
# which_model_infile = '../../ngc253/muse/2024July2/BIC_test/BIC_RATIO_1p5.fits'
which_model_infile = '/Users/serenac/Desktop/research/ngc253/muse/2024July2/BIC_RATIO_1_1p5_VelSep.fits'
hdu_BIC = fits.open(which_model_infile)
global which_model; which_model = hdu_BIC[0].data

# open the disk velocity model
disk_model_infile = '/Users/serenac/Desktop/research/ngc253/muse/data/ngc253_se_halpha_vel_model_smooth_replaceCO_andsmooth.fits'
#disk_model_infile = '../../ngc253/muse/data/NGC253_NKrieger_model.fits'
hdu_Ha = fits.open(disk_model_infile)
global disk_model; disk_model = hdu_Ha[0].data #- 243. ## FIXME: this is just for using the CO model