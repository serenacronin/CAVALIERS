import numpy as np
# import pandas as pd
import sys
sys.path.append('SeparateComps/')
import separateComps

def line_flux(amp, sig):
    
    """
    
    Adopted from https://lukeholden.com/blog/measuring_emission_line_fluxes_python.html#measuring_line_fluxes
    
    """
    
    sig_Ang = (sig * 2.2) / 100
    
    flux = amp * sig_Ang * np.sqrt(2*np.pi)
#     flux_err = flux*np.sqrt((peak_err/peak)**2 + (width_err/width)**2)
    return flux


def calcRatio(model):

    if model == 1:
        ha_disk_amp = separateComps.separateComps('amp', model, 'ha')
        nii_disk_amp = separateComps.separateComps('amp', model, 'niib')
        ha_disk_sig = separateComps.separateComps('sigma', model, 'ha')
        nii_disk_sig = separateComps.separateComps('sigma', model, 'niib')
        ratio = line_flux(nii_disk_amp, nii_disk_sig) / line_flux(ha_disk_amp, ha_disk_sig)
        return ratio

    if model == 2:
        ha_outflow_amp, ha_disk_amp = separateComps.separateComps('amp', model, 'ha')
        nii_outflow_amp, nii_disk_amp = separateComps.separateComps('amp', model, 'niib')
        ha_outflow_sig, ha_disk_sig = separateComps.separateComps('sigma', model, 'ha')
        nii_outflow_sig, nii_disk_sig = separateComps.separateComps('sigma', model, 'niib')
        ratio_disk = line_flux(nii_disk_amp, nii_disk_sig) / line_flux(ha_disk_amp, ha_disk_sig)
        ratio_outflow = line_flux(nii_outflow_amp, nii_outflow_sig) / line_flux(ha_outflow_amp, ha_outflow_sig)
        return ratio_outflow, ratio_disk
    
    if model == 3:
        ha_outflow1_amp, ha_outflow2_amp, ha_disk_amp = separateComps.separateComps('amp', model, 'ha')
        nii_outflow1_amp, nii_outflow2_amp, nii_disk_amp = separateComps.separateComps('amp', model, 'niib')
        ha_outflow1_sig, ha_outflow2_sig, ha_disk_sig = separateComps.separateComps('sigma', model, 'ha')
        nii_outflow1_sig, nii_outflow2_sig, nii_disk_sig = separateComps.separateComps('sigma', model, 'niib')
        ratio_disk = line_flux(nii_disk_amp, nii_disk_sig) / line_flux(ha_disk_amp, ha_disk_sig)
        ratio_outflow1 = line_flux(nii_outflow1_amp, nii_outflow1_sig) / line_flux(ha_outflow1_amp, ha_outflow1_sig)
        ratio_outflow2 = line_flux(nii_outflow2_amp, nii_outflow2_sig) / line_flux(ha_outflow2_amp, ha_outflow2_sig)
        return ratio_outflow1, ratio_outflow2, ratio_disk

    if model == 123:
        ha_disk_amp, ha_outflow1_amp, ha_outflow2_amp = separateComps.buildMaps('amp', line='ha')
        ha_disk_sig, ha_outflow1_sig, ha_outflow2_sig = separateComps.buildMaps('sigma', line='ha')
        nii_disk_amp, nii_outflow1_amp, nii_outflow2_amp = separateComps.buildMaps('amp', line='niib')
        nii_disk_sig, nii_outflow1_sig, nii_outflow2_sig = separateComps.buildMaps('sigma', line='niib')
        ratio_disk = line_flux(nii_disk_amp, nii_disk_sig) / line_flux(ha_disk_amp, ha_disk_sig)
        ratio_outflow1 = line_flux(nii_outflow1_amp, nii_outflow1_sig) / line_flux(ha_outflow1_amp, ha_outflow1_sig)
        ratio_outflow2 = line_flux(nii_outflow2_amp, nii_outflow2_sig) / line_flux(ha_outflow2_amp, ha_outflow2_sig)
        return ratio_outflow1, ratio_outflow2, ratio_disk
