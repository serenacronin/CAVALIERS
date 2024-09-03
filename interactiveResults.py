#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 9 16:18:50 2022
Updated Fri Aug 2 2024

@author: Serena A. Cronin

This script allows me to interact with my ratio plots!
I can now pull out a spectrum in each pixel with the click of a button!

"""

import sys
sys.path.append('../astro_tools')

from gauss_tools import one_gaussian, two_gaussian # type: ignore
import pandas as pd # type: ignore
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from astropy.io import fits # type: ignore
from matplotlib.lines import Line2D # type: ignore
from astropy.wcs import wcs # type: ignore
from spectral_cube import SpectralCube # type: ignore
from astropy import units as u # type: ignore
from matplotlib import interactive  # type: ignore
import cmasher as cmr # type:ignore
sys.path.append('SeparateComps/')
import getInfo

def build_gauss(x_data, line, model, fits):

    if model == 1:
        if line == 'niia':
            comp_num = 1
        elif line == 'ha':
            comp_num = 2
        elif line == 'niib':
            comp_num = 3

        comp = one_gaussian(x_data,
                    float(fits['Amp%s' % comp_num]),
                    float(fits['Wvl%s' % comp_num]),
                    float(fits['Sig%s' % comp_num]))
        return comp

    if model == 2:
        if line == 'niia':
            comp_num_blue = 1
            comp_num_red = 2
        elif line == 'ha':
            comp_num_blue = 3
            comp_num_red = 4
        elif line == 'niib':
            comp_num_blue = 5
            comp_num_red = 6

        comp_blue = one_gaussian(x_data, float(fits['Amp%s' % comp_num_blue]), 
                           float(fits['Wvl%s' % comp_num_blue]), 
                           float(fits['Sig%s' % comp_num_blue]))
    
        comp_red = one_gaussian(x_data,
                                float(fits['Amp%s' % comp_num_red]),
                                float(fits['Wvl%s' % comp_num_red]),
                                float(fits['Sig%s' % comp_num_red]))
        return comp_blue+comp_red, comp_blue, comp_red
    
    if model == 3:
        if line == 'niia':
            comp_num_blue = 1
            comp_num_mid = 2
            comp_num_red = 3
        elif line == 'ha':
            comp_num_blue = 4
            comp_num_mid = 5
            comp_num_red = 6
        elif line == 'niib':
            comp_num_blue = 7
            comp_num_mid = 8
            comp_num_red = 9

        comp_blue = one_gaussian(x_data, float(fits['Amp%s' % comp_num_blue]), 
                           float(fits['Wvl%s' % comp_num_blue]), 
                           float(fits['Sig%s' % comp_num_blue]))
        
        comp_mid = one_gaussian(x_data, float(fits['Amp%s' % comp_num_mid]), 
                           float(fits['Wvl%s' % comp_num_mid]), 
                           float(fits['Sig%s' % comp_num_mid]))
    
        comp_red = one_gaussian(x_data,
                                float(fits['Amp%s' % comp_num_red]),
                                float(fits['Wvl%s' % comp_num_red]),
                                float(fits['Sig%s' % comp_num_red]))
        return comp_blue+comp_mid+comp_red, comp_blue, comp_mid, comp_red


def plot_spec(x_data, spectrum, ix, iy):

    # make the plot!
    fig, ax = plt.subplots()
    ax.step(x_data, spectrum, color='gray', alpha=0.5, lw=2)
    
    # set axes labels
    ax.tick_params(axis='both', which='both',direction='in',
                   width=2.5, labelsize=16, length=7)
    ax.set_xlabel('Wavelength ($\AA$)', fontsize=20)
    ax.set_ylabel(r'Flux $(10^{-20} \mathrm{erg/s/cm^2/\AA})$', fontsize=20)
    
    # set a title
    ax.set_title('Pixel: (%s, %s)' % (int(round(ix,0)),int(round(iy,0))), fontsize=16)

    return plt

def plot_model1(x_data, spectrum, ix, iy):

    fits1 = pd.read_csv('../ngc253/muse/2024July2/fits1_total/fits1_reordered.txt')
    fits1_pix = fits1[(fits1['X'] == int(round(ix,0))) & (fits1['Y'] == int(round(iy,0)))]
    niib_tot = build_gauss(x_data, 'niib', 1, fits1_pix)
    ha_tot = build_gauss(x_data, 'ha', 1, fits1_pix)
    niia_tot = build_gauss(x_data, 'niia', 1, fits1_pix)

    # make the plot!
    fig, ax = plt.subplots()
    ax.step(x_data, spectrum, color='gray', alpha=0.5, lw=2)
    
    # plot the fits
    ax.plot(x_data, ha_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niib_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niia_tot, color='tab:pink', lw=2)
    
    # set axes labels
    ax.tick_params(axis='both', which='both',direction='in',
                   width=2.5, labelsize=16, length=7)
    ax.set_xlabel('Wavelength ($\AA$)', fontsize=20)
    ax.set_ylabel(r'Flux $(10^{-20} \mathrm{erg/s/cm^2/\AA})$', fontsize=20)
    
    # set a title
    ax.set_title('Pixel: (%s, %s)' % (int(round(ix,0)),int(round(iy,0))), fontsize=16)
    
    # plot a legend
    custom_lines = [Line2D([0], [0], color='tab:pink', lw=2),
                    Line2D([0], [0], color='white', lw=1),
                    Line2D([0], [0], color='white', lw=1)]
    plt.legend(custom_lines,['1 Component', 
                             'RedChiSq: %s' % np.round(float(fits1_pix['RedChiSq']),2), 
                             'BIC: %s' % np.round(float(fits1_pix['BIC']),2)],
                             fontsize=12, loc='upper left')

    return plt


def plot_model2(x_data, spectrum, ix, iy):

    fits2 = pd.read_csv('../ngc253/muse/2024July2/fits2_total/fits2_reordered.txt')
    fits2_pix = fits2[(fits2['X'] == int(round(ix,0))) & (fits2['Y'] == int(round(iy,0)))]
    niib_tot, niib_blue, niib_red = build_gauss(x_data, 'niib', 2, fits2_pix)
    ha_tot, ha_blue, ha_red = build_gauss(x_data, 'ha', 2, fits2_pix)
    niia_tot, niia_blue, niia_red = build_gauss(x_data, 'niia', 2, fits2_pix)

    # make the plot!
    fig, ax = plt.subplots()
    ax.step(x_data, spectrum, color='gray', alpha=0.5, lw=2)
    
    # plot the fits
    ax.plot(x_data, ha_tot, color='tab:pink', lw=2)
    ax.plot(x_data, ha_blue, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, ha_red, color = 'tab:cyan', lw=1)
    
    ax.plot(x_data, niib_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niib_blue, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, niib_red, color = 'tab:cyan', lw=1)

    ax.plot(x_data, niia_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niia_blue, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, niia_red, color = 'tab:cyan', lw=1)
    
    # set axes labels
    ax.tick_params(axis='both', which='both',direction='in',
                   width=2.5, labelsize=16, length=7)
    ax.set_xlabel('Wavelength ($\AA$)', fontsize=20)
    ax.set_ylabel(r'Flux $(10^{-20} \mathrm{erg/s/cm^2/\AA})$', fontsize=20)
    
    # set a title
    ax.set_title('Pixel: (%s, %s)' % (int(round(ix,0)),int(round(iy,0))), fontsize=16)
    
    # plot a legend
    custom_lines = [Line2D([0], [0], color='tab:pink', lw=2),
                    Line2D([0], [0], color='tab:cyan', lw=2),
                    Line2D([0], [0], color='white', lw=1),
                    Line2D([0], [0], color='white', lw=1)]
    plt.legend(custom_lines,['Composite', '2 Components', 
                             'RedChiSq: %s' % np.round(float(fits2_pix['RedChiSq']),2), 
                             'BIC: %s' % np.round(float(fits2_pix['BIC']),2)],
                             fontsize=12, loc='upper left')

    return plt


def plot_model3(x_data, spectrum, ix, iy):

    fits3 = pd.read_csv('../ngc253/muse/2024July2/fits3_total/fits3_reordered.txt')
    fits3_pix = fits3[(fits3['X'] == int(round(ix,0))) & (fits3['Y'] == int(round(iy,0)))]
    niib_tot, niib_blue, niib_mid, niib_red = build_gauss(x_data, 'niib', 3, fits3_pix)
    ha_tot, ha_blue, ha_mid, ha_red = build_gauss(x_data, 'ha', 3, fits3_pix)
    niia_tot, niia_blue, niia_mid, niia_red = build_gauss(x_data, 'niia', 3, fits3_pix)

    # make the plot!
    fig, ax = plt.subplots()
    ax.step(x_data, spectrum, color='gray', alpha=0.5, lw=2)
    
    # plot the fits
    ax.plot(x_data, ha_tot, color='tab:pink', lw=2)
    ax.plot(x_data, ha_blue, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, ha_mid, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, ha_red, color = 'tab:cyan', lw=1)
    
    ax.plot(x_data, niib_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niib_blue, color = 'tab:cyan', lw=1)
    ax.plot(x_data, niib_mid, color = 'tab:cyan', lw=1)  
    ax.plot(x_data, niib_red, color = 'tab:cyan', lw=1)

    ax.plot(x_data, niia_tot, color='tab:pink', lw=2)
    ax.plot(x_data, niia_blue, color = 'tab:cyan', lw=1)
    ax.plot(x_data, niia_mid, color = 'tab:cyan', lw=1) 
    ax.plot(x_data, niia_red, color = 'tab:cyan', lw=1)
    
    # set axes labels
    ax.tick_params(axis='both', which='both',direction='in',
                   width=2.5, labelsize=16, length=7)
    ax.set_xlabel('Wavelength ($\AA$)', fontsize=20)
    ax.set_ylabel(r'Flux $(10^{-20} \mathrm{erg/s/cm^2/\AA})$', fontsize=20)
    
    # set a title
    ax.set_title('Pixel: (%s, %s)' % (int(round(ix,0)),int(round(iy,0))), fontsize=16)
    
    # plot a legend
    custom_lines = [Line2D([0], [0], color='tab:pink', lw=2),
                    Line2D([0], [0], color='tab:cyan', lw=2),
                    Line2D([0], [0], color='white', lw=1),
                    Line2D([0], [0], color='white', lw=1)]
    ax.legend(custom_lines,['Composite', '3 Components', 
                             'RedChiSq: %s' % np.round(float(fits3_pix['RedChiSq']),2), 
                             'BIC: %s' % np.round(float(fits3_pix['BIC']),2)],
                             fontsize=12, loc='upper left')

    return plt


def onpick(event):
    
    """
    
    Function that determines what to do when the user clicks
    the plot!
    
    event : Class
        Stores the event information (e.g., a mouse click)
    ParCubeFile : str
        File string of the data cube that contains the fit parameter info
        needed for the pop-up plots.
    
    """
    
    # spectral cube for extracting the spectrum per pixel
    # variables for my CreateCube function

    SlabLower = 6500
    SlabUpper = 6650
    filename = '../ngc253/muse/data/NGC253_MUSE_SE_Fitted.fits'
    SpecCube = SpectralCube.read(filename, hdu=1).spectral_slab(SlabLower * u.AA, 
																SlabUpper * u.AA)
      
    # THE ORDER OF IX, IY IS CORRECT FOR THE NEXT SEVERAL LINES
    # IDK WHY IT IS DOING THIS I GUESS EVENT SWAPPED THE AXES
    # IDK
    # JUST DON'T CHANGE THEM A;LSDKFKL;JASDL;KFJKSLA
    # get the coordinates of the event
    ix, iy = event.xdata, event.ydata
    
    # plot the actual spectrum from the spectral cube
    # get the x and y-axis data
    spectrum = np.array(SpecCube[:,int(round(iy,0)),int(round(ix,0))], dtype='float64')  # y-axis data
    minval = min(np.array(SpecCube.spectral_axis))
    maxval = max(np.array(SpecCube.spectral_axis))
    x_data = np.linspace(minval, maxval, len(spectrum)) # x-axis data

    plt = plot_spec(x_data, spectrum, ix, iy)
    interactive(True)
    plt.show()

    plt = plot_model1(x_data, spectrum, ix, iy)
    interactive(True)
    plt.show()

    plt = plot_model2(x_data, spectrum, ix, iy)
    interactive(True)
    plt.show()

    plt = plot_model3(x_data, spectrum, ix, iy)
    interactive(False)
    plt.show()

    return

##### --- SET UP PLOT --- #####
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2.5
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.variant"] = "small-caps"

BIC = '../ngc253/muse/2024July2/BIC_RATIO_1p25_1p5_VelSep.fits'
vmin = 1
vmax = 3
cmap = cmr.get_sub_cmap('cmr.cosmic', 0.2, 1, N=3)
label = 'Num. of Components'
title = 'Click me!'

# fits data cube file
BIC_data = fits.open(BIC)[0]

# plot the map
fig, ax = plt.subplots(1,1)
im = ax.imshow(BIC_data.data, vmin=vmin, vmax=vmax, 
                origin='lower', cmap=cmap, picker=True)

ax.tick_params(axis='both', which='both',direction='in',
               width=2.5, labelsize=16, length=7)
ax.set_title(title, fontsize=16)

# ax.set_xlabel('R.A. Offset from %s (arcsec)' % getInfo.offset_ra, fontsize=14)
# ax.set_ylabel('Dec. Offset from %s (arcsec)' % getInfo.offset_dec, fontsize=14)
 
bar = plt.colorbar(im, fraction=0.046)
bar.set_label(label, fontsize=18)
bar.ax.tick_params(width=2.5, labelsize=16, length=7, direction='in')
bar.ax.locator_params(nbins=3)

fig.canvas.mpl_connect('button_press_event', onpick)
plt.show()