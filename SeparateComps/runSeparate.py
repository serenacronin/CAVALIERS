import sys
sys.path.append('SeparateComps/')
import separateComps, plotComps, getInfo, lineRatios
from matplotlib import interactive  # type: ignore

param = 'vel'
model = 123
line = 'niib'
cbar_range = [-300, 300]
filt_bad_vals = False
save_as_fits = True

if model == 1:
    disk_map = separateComps.separateComps(param, model, line)
    plt = plotComps.plotComps(disk_map, comp='Disk', param=param, model=1, 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)  # allows both plots to show up
    plt.show()

    plt = plotComps.plotWhichModel(extent=getInfo.extent,
                                   offset_ra=getInfo.offset_ra,
                                   offset_dec=getInfo.offset_dec)
    interactive(True)  # allows both plots to show up
    plt.show()

    # ratio
    ratio = lineRatios.calcRatio(model)
    plt = plotComps.plotRatio(ratio, 'Disk', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(False)
    plt.show()



elif model == 2:
    outflow_map, disk_map = separateComps.separateComps(param, model, line)

    plt = plotComps.plotComps(outflow_map, comp='Outflow', param=param, model=2, 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotComps(disk_map, comp='Disk', param=param, model=2, 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotWhichModel(extent=getInfo.extent, 
                                   offset_ra=getInfo.offset_ra,
                                   offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    # ratio
    ratio_outflow, ratio_disk = lineRatios.calcRatio(model)
    plt = plotComps.plotRatio(ratio_outflow, 'Outflow', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotRatio(ratio_disk, 'Disk', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(False)
    plt.show()


elif model == 3:
    outflow_map1, outflow_map2, disk_map = separateComps.separateComps(param, model, line)

    plt = plotComps.plotComps(outflow_map1, comp='Outflow1', param=param, model=3, 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotComps(outflow_map2, comp='Outflow2', param=param, model=3, 
                    cbar_range=cbar_range, extent=getInfo.extent,
                    offset_ra=getInfo.offset_ra, 
                    offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotComps(disk_map, comp='Disk', param=param, model=3, 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotWhichModel(extent=getInfo.extent, 
                                   offset_ra=getInfo.offset_ra,
                                   offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    # ratio
    ratio_outflow1, ratio_outflow2, ratio_disk = lineRatios.calcRatio(model)
    plt = plotComps.plotRatio(ratio_outflow1, 'Outflow1', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotRatio(ratio_outflow2, 'Outflow2', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotRatio(ratio_disk, 'Disk', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(False)
    plt.show()


elif model == 123:

    if filt_bad_vals == True:
        disk_map_total, outflow_map1_total, outflow_map2_total = separateComps.filtered_maps(param, line)
    else:
        disk_map_total, outflow_map1_total, outflow_map2_total = separateComps.buildMaps(param, line=line)

    plt = plotComps.plotComps(outflow_map1_total, comp='Outflow1', param=param, model='all', 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotComps(outflow_map2_total, comp='Outflow2', param=param, model='all', 
                    cbar_range=cbar_range, extent=getInfo.extent,
                    offset_ra=getInfo.offset_ra, 
                    offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotComps(disk_map_total, comp='Disk', param=param, model='all', 
                        cbar_range=cbar_range, extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra, 
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotWhichModel(extent=getInfo.extent, 
                                   offset_ra=getInfo.offset_ra,
                                   offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    # ratio
    ratio_outflow1, ratio_outflow2, ratio_disk = lineRatios.calcRatio(model)
    plt = plotComps.plotRatio(ratio_outflow1, 'Outflow1', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotRatio(ratio_outflow2, 'Outflow2', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(True)
    plt.show()

    plt = plotComps.plotRatio(ratio_disk, 'Disk', model, [0,2], extent=getInfo.extent,
                        offset_ra=getInfo.offset_ra,
                        offset_dec=getInfo.offset_dec)
    interactive(False)
    plt.show()

    # save everything as a fits file if we are satisfied
    if save_as_fits == True:
        disk_map_total, outflow_map1_total, outflow_map2_total
        ratio_disk, ratio_outflow1, ratio_outflow2
