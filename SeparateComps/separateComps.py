import numpy as np
import sys
sys.path.append('SeparateComps/')
import getInfo

# access global variables from getInfo.py
global fits1; fits1 = getInfo.fits1
global fits2; fits2 = getInfo.fits2
global fits3; fits3 = getInfo.fits3
global y; y = getInfo.y
global x; x = getInfo.x
global which_model; which_model = getInfo.which_model
global disk_model; disk_model = getInfo.disk_model

def separateComps(param, model, line):

    """
    comp : which component do you want to look at? 
            Outflow Blue, Outflow Red, Disk
    param : which parameter do you want to look at?
            vel, amp, sigma
    model : which model (i.e., number of lines)?
            1, 2, 3
    
    """
    
    ## One component fit, where the map will just be the disk map ##
    if model == 1:
        if line == 'niib':
            comp_num = 3
        elif line == 'ha':
            comp_num = 2
        param_map = np.full((y,x), np.nan)
        for index, row in fits1.iterrows():
            amp = row['Amp%s' % comp_num]
            if float(amp) < 5.:  # if the component amp is close to 0.
                param_map[int(row['Y']), int(row['X'])] = np.nan  # blank out that component

            else:
                if param == 'BIC':
                    param_map[int(row['Y']), int(row['X'])] = row['BIC']
                if param == 'amp':
                    param_map[int(row['Y']), int(row['X'])] = amp
                if param == 'vel':
                    param_map[int(row['Y']), int(row['X'])] = row['Vel%s' % comp_num]
                if param == 'sigma':
                    param_map[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num]
        
        disk_map = np.full((y,x), fill_value = np.nan)
        disk_map[which_model == 1.] = param_map[which_model == 1.]

        return disk_map

    ## Two component fit, where the maps will be disk map and outflow map ##
    if model == 2:
        if line == 'niib':
            comp_num_blue = 5
            comp_num_red = 6
        elif line == 'ha':
            comp_num_blue = 3
            comp_num_red = 4
        vel_map_blue = np.full((y,x), np.nan)  # need for splitting on vels
        vel_map_red = np.full((y,x), np.nan) 
        param_map_blue = np.full((y,x), np.nan)  # if the parameter is not vel
        param_map_red = np.full((y,x), np.nan)

        for index, row in fits2.iterrows():

            amp_blue = row['Amp%s' % comp_num_blue]
            amp_red = row['Amp%s' % comp_num_red]
            if float(amp_blue) < 5.:  # if the component amp is close to 0.
                param_map_blue[int(row['Y']), int(row['X'])] = 1e9  # blank out that component
            elif float(amp_red) < 5.:  # if the component amp is close to 0.
                param_map_red[int(row['Y']), int(row['X'])] = 1e9  # blank out that component

            else:
                # grab the velocities to separate components based on them
                vel_map_blue[int(row['Y']), int(row['X'])] = row['Vel5']
                vel_map_red[int(row['Y']), int(row['X'])] = row['Vel6']

                # grab the parameters you want to plot
                if param == 'BIC':
                    param_map_blue[int(row['Y']), int(row['X'])] = row['BIC']
                    param_map_red[int(row['Y']), int(row['X'])] = row['BIC']
                if param == 'amp':
                    param_map_blue[int(row['Y']), int(row['X'])] = row['Amp%s' % comp_num_blue]
                    param_map_red[int(row['Y']), int(row['X'])] = row['Amp%s' % comp_num_red]
                if param == 'sigma':
                    param_map_blue[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num_blue]
                    param_map_red[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num_red]
                if param == 'vel':
                    param_map_blue = vel_map_blue
                    param_map_red = vel_map_red

        # blueshifted 2 component model is disk
        mask_blue = ((np.abs(vel_map_blue - disk_model) < np.abs(vel_map_red - disk_model)))
        # redshifted 2 component model is disk
        mask_red = ((np.abs(vel_map_red - disk_model) < np.abs(vel_map_blue - disk_model)))

        # apply the mask to separate the components
        disk_map = np.full((y,x), fill_value = np.nan)
        disk_map[(which_model == 2.) & mask_blue] = \
            param_map_blue[(which_model == 2.) & mask_blue]
        disk_map[(which_model == 2.) & mask_red] = \
            param_map_red[(which_model == 2.) & mask_red]

        outflow_map = np.full((y,x), fill_value = np.nan)
        outflow_map[(which_model == 2.) & mask_red] = \
            param_map_blue[(which_model == 2.) & mask_red]
        outflow_map[(which_model == 2.) & mask_blue] = \
            param_map_red[(which_model == 2.) & mask_blue]
        
        outflow_map[outflow_map == 1e6] = np.nan
        disk_map[disk_map == 1e6] = np.nan
        
        return outflow_map, disk_map
    

    ## Three component fit, where the maps will be disk map,
    ##  blue outflow map, and red outflow map ##
    if model == 3:
        if line == 'niib':
            comp_num_blue = 7
            comp_num_mid = 8
            comp_num_red = 9
        elif line == 'ha':
            comp_num_blue = 4
            comp_num_mid = 5
            comp_num_red = 6
        vel_map_blue = np.full((y,x), np.nan)  # need for splitting on vels
        vel_map_mid = np.full((y,x), np.nan) 
        vel_map_red = np.full((y,x), np.nan) 
        param_map_blue = np.full((y,x), np.nan)  # if the parameter is not vel
        param_map_mid = np.full((y,x), np.nan)
        param_map_red = np.full((y,x), np.nan)

        for index, row in fits3.iterrows():

            # grab the velocities to separate components based on them
            vel_map_blue[int(row['Y']), int(row['X'])] = row['Vel%s' % comp_num_blue]
            vel_map_mid[int(row['Y']), int(row['X'])] = row['Vel%s' % comp_num_mid]
            vel_map_red[int(row['Y']), int(row['X'])] = row['Vel%s' % comp_num_red]

            # grab the parameters you want to plot
            if param == 'BIC':
                param_map_blue[int(row['Y']), int(row['X'])] = row['BIC']
                param_map_red[int(row['Y']), int(row['X'])] = row['BIC']
            if param == 'amp':
                param_map_blue[int(row['Y']), int(row['X'])] = row['Amp%s' % comp_num_blue]
                param_map_mid[int(row['Y']), int(row['X'])] = row['Amp%s' % comp_num_mid]
                param_map_red[int(row['Y']), int(row['X'])] = row['Amp%s' % comp_num_red]
            if param == 'sigma':
                param_map_blue[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num_blue]
                param_map_mid[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num_mid]
                param_map_red[int(row['Y']), int(row['X'])] = row['SigVel%s' % comp_num_red]

        # if the parameter to plot is velocity, we already have it!
        if param == 'vel':
            param_map_blue = vel_map_blue
            param_map_mid = vel_map_mid
            param_map_red = vel_map_red

        # blueshifted 3 component model is disk
        mask_blue = ((np.abs(vel_map_blue - disk_model) < np.abs(vel_map_red - disk_model)) & \
                   (np.abs(vel_map_blue - disk_model) < np.abs(vel_map_mid - disk_model)))

        # middle 3 component model is disk
        mask_mid = ((np.abs(vel_map_mid - disk_model) < np.abs(vel_map_blue - disk_model)) & \
                (np.abs(vel_map_mid - disk_model) < np.abs(vel_map_red - disk_model)))

        # redshifted 3 component model is disk
        mask_red = ((np.abs(vel_map_red - disk_model) < np.abs(vel_map_blue - disk_model)) & \
                (np.abs(vel_map_red - disk_model) < np.abs(vel_map_mid - disk_model)))

        # apply the mask to separate the components
        disk_map = np.full((y,x), fill_value = np.nan)
        disk_map[(which_model == 3.) & mask_blue] = \
            param_map_blue[(which_model == 3.) & mask_blue]
        disk_map[(which_model == 3.) & mask_mid] = \
            param_map_mid[(which_model == 3.) & mask_mid]
        disk_map[(which_model == 3.) & mask_red] = \
            param_map_red[(which_model == 3.) & mask_red]
        
        # if the disk is redshifted or the middle component
        # then the outflow is blueshifted; if the disk is
        # blueshifted, then the disk is the middle component
        outflow_map1 = np.full((y,x), fill_value = np.nan)
        outflow_map1[(which_model == 3.) & mask_red] = \
            param_map_blue[(which_model == 3.) & mask_red]
        outflow_map1[(which_model == 3.) & mask_mid] = \
            param_map_blue[(which_model == 3.) & mask_mid]
        outflow_map1[(which_model == 3.) & mask_blue] = \
            param_map_mid[(which_model == 3.) & mask_blue]
        
        # and the opposite here
        outflow_map2 = np.full((y,x), fill_value = np.nan)
        outflow_map2[(which_model == 3.) & mask_red] = \
            param_map_mid[(which_model == 3.) & mask_red]
        outflow_map2[(which_model == 3.) & mask_blue] = \
            param_map_red[(which_model == 3.) & mask_blue]
        outflow_map2[(which_model == 3.) & mask_mid] = \
            param_map_red[(which_model == 3.) & mask_mid]
        
        return outflow_map1, outflow_map2, disk_map

def buildMaps(param, line):
    disk_map_model1 = separateComps(param, model=1, line=line)
    outflow_map_model2, disk_map_model2 = separateComps(param, model=2, line=line)
    outflow_map1_model3, outflow_map2_model3, disk_map_model3 = separateComps(param, model=3, line=line)
    
    disk_map_total = np.full(disk_map_model1.shape, fill_value=np.nan)
    outflow_map1_total = np.full(disk_map_model1.shape, fill_value=np.nan)
    outflow_map2_total = np.full(disk_map_model1.shape, fill_value=np.nan)

    disk_map_model1[~np.isfinite(disk_map_model1)] = 0.
    disk_map_model2[~np.isfinite(disk_map_model2)] = 0.
    disk_map_model3[~np.isfinite(disk_map_model3)] = 0.
    disk_map_total = disk_map_model1 + disk_map_model2 + disk_map_model3
    disk_map_total[outflow_map1_total == 0.] = np.nan

    outflow_map_model2[~np.isfinite(outflow_map_model2)] = 0.
    outflow_map1_model3[~np.isfinite(outflow_map1_model3)] = 0.
    outflow_map1_total = outflow_map_model2 + outflow_map1_model3
    outflow_map1_total[outflow_map1_total == 0.] = np.nan

    outflow_map2_total = outflow_map2_model3

    return(disk_map_total, outflow_map1_total, outflow_map2_total)


def interp_bad_values(mapp, vel_mapp, thresh):

    # find indices with values above the vel threshold
    bad_indices = np.argwhere(np.abs(vel_mapp) > thresh)
    filtered_mapp = mapp.copy()

    for i in range(len(bad_indices)):
        # good_y = (bad_indices[i][1]+1, bad_indices[i][1]+1, bad_indices[i][1]-1, bad_indices[i][1]-1)
        # good_x = (bad_indices[i][0]+1, bad_indices[i][0]+1, bad_indices[i][0]-1, bad_indices[i][0]-1)
        filt = mapp[(bad_indices[i][0]-5):(bad_indices[i][0]+5), (bad_indices[i][1]-5):(bad_indices[i][1]+5)]
        filt_vel = vel_mapp[(bad_indices[i][0]-5):(bad_indices[i][0]+5), (bad_indices[i][1]-5):(bad_indices[i][1]+5)]
        good_val = np.mean(filt[np.abs(filt_vel) < thresh])
        filtered_mapp[bad_indices[i][0], bad_indices[i][1]] = good_val

    return(filtered_mapp)


def filtered_maps(param, line):

    # need to filter based on velocity outliers
    disk_vel_map, outflow1_vel_map, outflow2_vel_map = buildMaps('vel', line)

    # if we only want velocity maps
    if param == 'vel':
        disk_map_total = disk_vel_map.copy()
        outflow1_map_total = outflow1_vel_map.copy()
        outflow2_map_total = outflow2_vel_map.copy()
    # if we want another parameter, like amp, build those maps
    else:
        disk_map_total, outflow1_map_total, outflow2_map_total = buildMaps(param, line)
    
    # filter the values
    disk_map_filtered = interp_bad_values(disk_map_total, disk_vel_map, thresh=900)
    outflow1_map_filtered = interp_bad_values(outflow1_map_total, outflow1_vel_map, thresh=900)
    outflow2_map_filtered = interp_bad_values(outflow2_map_total, outflow2_vel_map, thresh=900)

    return(disk_map_filtered, outflow1_map_filtered, outflow2_map_filtered)