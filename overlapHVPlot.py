#!/usr/bin/env python3
'''
File to Graph Overlapping Plot. created as a seperate file due to key differences in the way it runs from mkHVScanPlot.
'''
import numpy as np
import matplotlib.pyplot as plt
import src.DataFile as df
import src.helpers as hp

def mkOverlappingPlot(mscan_list1, strip1, src1, hole1, mscan_list2, **kwargs):       
    '''
        Plot overlapping Currents given two mscan_lists. Currently expects two files of same strip, src, and hole   
    '''
    defaults = {
        '''
            Holds more parameters. More documentation in mkPlateauPlot if needed.
        '''
        #optional parameters
        'start_point': 0,           #(V), starting point in data run.
        'exclude_start' : True,     #exludes starting point of data(intended to remove zero)
        'end_point': 500,           #(V), ending point in data
        #change scale of graph image
        'x_low_lim': None,          
        'x_upper_lim': None,
        'y_lower_lim': None,
        'y_upper_lim': None,

    }
    #set parameters in
    defaults.update(kwargs)
    
    #Start and stop graphing points for first data file
    first_start_fit = np.where(mscan_list1[0] == 0)[0][0]
    first_stop_fit = np.where(mscan_list1[0] == 500)[0][0]+1

    #start and stop graphing points for second data file
    second_start_fit = np.where(mscan_list2[0] == 0)[0][0]
    second_stop_fit = np.where(mscan_list2[0] == 500)[0][0]+1

    #Whether to exclude first point or not(0V point.) defaults to true unless exclude_start parameter is given
    first_start_fit += 1 if defaults['exclude_start'] == True else None
    second_start_fit += 1 if defaults['exclude_start'] == True else None

    #get x and y vals for first data file
    first_x_vals = mscan_list1[0][first_start_fit:first_stop_fit]
    first_y_vals = mscan_list1[1][first_start_fit:first_stop_fit]

    #get x and y vals for second data file
    second_x_vals = mscan_list2[0][second_start_fit:second_stop_fit]
    second_y_vals = mscan_list2[1][second_start_fit:second_stop_fit]


    #plot both sets of errorbars 
    plt.errorbar(first_x_vals, first_y_vals, marker='.', yerr=mscan_list1[2][first_start_fit:first_stop_fit], linestyle='', label='first Avg Current')
    plt.errorbar(second_x_vals,second_y_vals, marker='.', yerr=mscan_list2[2][second_start_fit:second_stop_fit],linestyle='', label='Second Avg Current')

    #Calculate slope and intercept for both data files
    #slope1, intercept1 = np.polyfit(first_x_vals, first_y_vals, 1)
    #slope2, intercept2 = np.polyfit(second_x_vals, second_y_vals, 1)

    #Calculate both fit lines
    #first_fit_line = slope1 * first_x_vals + intercept1
    #second_fit_line = slope2 * second_x_vals + intercept2

    #Scatter points for both files
    #plt.scatter(first_x_vals, first_y_vals)
    #plt.scatter(second_x_vals, second_y_vals)

    #plot Linear fits for both files
    #plt.plot(first_x_vals, first_fit_line, color='red')
    #plt.plot(second_x_vals, second_fit_line, color='blue') 

    #Label Graph
    plt.title(f'MiniCSC4: HV Scan, L1, 90{src1}-Src1, {hole1.replace("_", "").replace("0","").upper()}, {strip1}, Run Comparison')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I')
    plt.legend()
    
    #Scaling
    plt.xlim(defaults['x_low_lim'], defaults['x_upper_lim'])
    plt.ylim(defaults['y_lower_lim'], defaults['y_upper_lim'])

    #show graph(for testing)
    plt.show()

    #Save and close graph
    #plt.savefig(f'./plots/HV_Scans/Run_Comparison/{strip1}_src{src1}_comparisonLow.png', format='png', dpi=400)
    plt.close()

def main(run_nm1, run_nm2):
    '''
        Main function for the file. Gets dataFile objs from createRun, matches data points, subtracts backgrounds, and makes 
    '''
    criteria = {'strip':int(run_nm1.split('S')[-1])} 
    strip = run_nm1.split("_")[-1]

    # Load HV Scans with source
    print(f'Loading Data for {run_nm1}')
    sdf1 = hp.createRun(run_nm1, criteria)
    sdf2 = hp.createRun(run_nm2, criteria)

     # Load Dark Current runs
    print(f'\nLoading Dark Currents for {run_nm1}')
    ddf1 = hp.createRun(run_nm1, criteria, src=False)
    ddf2 = hp.createRun(run_nm2, criteria, src=False)


     # Match the data points from the two scans
    print('\n\tMatching the HV points of the first scan')
    src_hvscan1 = sdf1.getHVScan()
    drk_hvscan1 = ddf1.getHVScan(src=False)
    src_hvscan1, drk_hvscan1 = hp.matching(src_hvscan1,drk_hvscan1)
    print('\n\tSubtracting background')
    corrected_hvscan1 = src_hvscan1[1] - drk_hvscan1[1]

    print("\n\tMatching second scan")
    src_hvscan2 = sdf2.getHVScan()
    drk_hvscan2 = ddf2.getHVScan(src=False)
    src_hvscan2, drk_hvscan2 = hp.matching(src_hvscan2, drk_hvscan2)
    print('\n\tSubtracting background 2nd run')
    corrected_hvscan2 = src_hvscan2[1] - drk_hvscan2[1]

    #NOTE: MORE NOISE THAN ERROR FUNCTIONALITY HAS BEEN REMOVED.


    #Holds HV, corrected Avg Current, combined stderror
    mhvscan1 = [src_hvscan1[0], corrected_hvscan1, hp.quadSum(src_hvscan1[2], drk_hvscan1[2])]
    src1 = sdf1.getFileSrc()
    hole1 = sdf1.getHole()

    mhvscan2 = [src_hvscan2[0], corrected_hvscan2, hp.quadSum(src_hvscan2[2], drk_hvscan2[2])]
    

    mkOverlappingPlot(mhvscan1, strip, src1, hole1, mhvscan2)

stripWdth = 12.7 #mm
stripGap = 0.35 ##mm <- ME1/1 chamber TDR (ME2/1 is 0.5 mm)

#NOTE: lists of different run files.comment and uncomment as needed. should be in a .data/HV_Scans folder.
#run_nms = ['241001_refMeasures_S2','241002_refMeasures_S6','240820_refMeasures_S7','241002_refMeasures_S10']
run_nm1 = '241001_refMeasures_S2'
run_nm2 = '241014_Plateau_S2'

main(run_nm1, run_nm2)



