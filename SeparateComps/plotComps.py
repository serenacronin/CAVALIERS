import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import cmasher as cmr
import sys
sys.path.append('SeparateComps/')
import getInfo
import numpy as np

def coneModel(x_apex=310, y_apex=320):

    """
    if you need to understand this 
    please see Emma Mirizio
    she won't know either
    she will have forgotten in 3 days
    and she apologizes for redefining variables frequently
    
    """

    m = np.tan(200*np.pi/180)
    b = y_apex - m*x_apex
    x1 = x_apex
    y1 = m*x1+b
    x2 = 40
    y2 = m*x2+b
    l1 = 15*np.cos(20*np.pi/180)/0.2
    x1 = x1 - l1
    y1 = m*x1+b
    solid_line1 = [[x1,x2], [y1,y2]]

    l2 = 35*np.cos(20*np.pi/180)/0.2
    x2 = x_apex - l2
    y2 = m*x2+b
    dotted_line1 = [[x1,x2], [y1,y2]]

    m = np.tan(260*np.pi/180)
    b = y_apex - m*x_apex
    x1 = x_apex
    y1 = m*x1+b
    x2 = 260
    y2 = m*x2+b
    l1 = 15*np.cos(80*np.pi/180)/0.2
    x1 = x1 - l1
    y1 = m*x1+b
    solid_line2 = [[x1,x2], [y1,y2]]

    l2 = 35*np.cos(80*np.pi/180)/0.2
    x2 = x_apex - l2
    y2 = m*x2+b
    dotted_line2 = [[x1,x2], [y1,y2]]

    return(solid_line1, dotted_line1, solid_line2, dotted_line2)



def plotComps(dat, comp, param, model, cbar_range, extent, offset_ra, offset_dec):

    """
    dat : the data map to plot
    comp : which component are we looking at? Outflow Blue, Red, or Disk?
    param : which parameter are we looking at? velocity, amplitude, or sigma?
    model : which fit did we do? 1, 2, or 3 components?
    plot_range : what range to plot
    
    """

    fig, ax = plt.subplots(1,1, figsize=(8,6))
    if param == 'vel':
        cmap = cmr.get_sub_cmap('cmr.fusion_r', 0.1, 0.90)

    im = ax.imshow(dat, origin='lower', cmap=cmap, vmin=cbar_range[0], vmax=cbar_range[1])
    ax.set_title('%s (%s-component model)' % (comp, model), fontsize=24)
    ax.tick_params(axis='both', which='both',direction='in', width=1.5, labelsize=18, length=5, color='white')
    ax.set_facecolor('black')

    at = AnchoredText(r'[N II]$\lambda6583$',  frameon=False, loc='lower left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    at = AnchoredText('SE', frameon=False, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    solid_line1, dotted_line1, solid_line2, dotted_line2 = coneModel()
    ax.plot(solid_line1[0], color='gray', lw=3)
    ax.plot(dotted_line1[0], color='gray', ls='--', lw=3)
    ax.plot(solid_line2[0], color='gray', lw=3)
    ax.plot(dotted_line2[0], color='gray', ls='--', lw=3)

    # ax.set_xlabel('R.A. Offset from %s (arcsec)' % offset_ra, fontsize=20)
    # ax.set_ylabel('Dec. Offset from %s (arcsec)' % offset_dec, fontsize=20)

    cax = ax.inset_axes([1.04, 0.002, 0.05, 1])  # make another axis for the colorbar 
                                                # [x0, y0, width, height] where x0, y0 = lower left corner
    bar = fig.colorbar(im, ax=ax, cax=cax)
    bar.set_label('%s' % param, fontsize=18)
    bar.ax.tick_params(width=2.5, labelsize=16, length=7, direction='in')

    plt.savefig(getInfo.savepath + 'se_%scomp_%s_%s.png' % (model, param, comp), dpi=200)
    return plt


def plotWhichModel(extent, offset_ra, offset_dec):

    """
    
    """

    which_model = getInfo.which_model

    # FIXME: way to add which exact BIC model I'm using in 
    # title and the file name

    fig, ax = plt.subplots(1,1, figsize=(8,6))
    cmap = cmr.get_sub_cmap('cmr.cosmic', 0.2, 1, N=3)
    im = ax.imshow(which_model, vmin=1, vmax=3, cmap=cmap, origin='lower')
    ax.tick_params(axis='both', which='both',direction='in', width=1.5, labelsize=18, length=5, color='white')

    # ax.set_xlabel('R.A. Offset from %s (arcsec)' % offset_ra, fontsize=20)
    # ax.set_ylabel('Dec. Offset from %s (arcsec)' % offset_dec, fontsize=20)

    cax = ax.inset_axes([1.04, 0.002, 0.05, 1])  # make another axis for the colorbar 
                                                # [x0, y0, width, height] where x0, y0 = lower left corner
    bar = fig.colorbar(im, ax=ax, cax=cax)
    bar.set_label('Number of Components', fontsize=18)
    bar.ax.tick_params(width=2.5, labelsize=16, length=7, direction='in')
    bar.ax.locator_params(nbins=3)

    plt.savefig(getInfo.savepath + 'which_model.png', dpi=200)
    return plt


def plotRatio(ratio, comp, model, cbar_range, extent, offset_ra, offset_dec):

    fig, ax = plt.subplots(1,1, figsize=(8,6))
    cmap = cmr.get_sub_cmap('cmr.tropical_r', 0.1, 0.90)
    im = ax.imshow(ratio, origin='lower', cmap=cmap, vmin=cbar_range[0], vmax=cbar_range[1])
    ax.set_title('%s (%s-component model)' % (comp, model), fontsize=24)
    ax.tick_params(axis='both', which='both',direction='in', width=1.5, labelsize=18, length=5, color='white')
    ax.set_facecolor('black')

    at = AnchoredText('SE', frameon=False, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    # ax.set_xlabel('R.A. Offset from %s (arcsec)' % offset_ra, fontsize=20)
    # ax.set_ylabel('Dec. Offset from %s (arcsec)' % offset_dec, fontsize=20)

    cax = ax.inset_axes([1.04, 0.002, 0.05, 1])  # make another axis for the colorbar 
                                                # [x0, y0, width, height] where x0, y0 = lower left corner
    bar = fig.colorbar(im, ax=ax, cax=cax)
    bar.set_label(r'[N II]$\lambda6583$ / H$\alpha$', fontsize=18)
    bar.ax.tick_params(width=2.5, labelsize=16, length=7, direction='in')
    bar.ax.axhline(10**(-0.2), c='w')

    plt.savefig(getInfo.savepath + 'se_%scomp_%s_nii_ha.png' % (model, comp), dpi=200)

    return plt