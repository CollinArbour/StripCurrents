import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.optimize import curve_fit

import helpers.DataFile as df
import helpers.functions as func
import helpers.plotting as mplt

def mkRawScans(strips,ax,mlabel=''):
    '''
    args:strips - [xpos,values,stderrs] **values is nA/mm (across with width of the strip)
    args:ax     -  plt axis object to plot on
    '''

    if mlabel:
        ax.errorbar(strips[0],strips[1],yerr=strips[2],linestyle='',marker='.',label=mlabel)
    else:
        ax.errorbar(strips[0],strips[1],yerr=strips[2],linestyle='',marker='.')
    
    ax.set_xlabel('Strip Position (mm)')
    ax.set_ylabel('Strip Current (nA)')

    return ax

def createRun(run_nm, criteria, src=True):
    '''Does Basic Data processing using the dataFile class. made into a function for code simplification.'''
    if src == False:
        dataObj = df.DataFile(f'{run_nm[7:]}_dark')
        dataObj.parseDataFileText(f'./data/HV_Scans/{run_nm}_dark.txt')
    else:
        dataObj = df.DataFile(f'{run_nm[7:]}')
        dataObj.parseDataFileText(f'./data/HV_Scans/{run_nm}.txt')
    
    dataObj.filterRuns(criteria)
    dataObj.sortDataRuns('hv')

    return dataObj

def matching(src,drk):
    '''
        This function takes in a source scan and dark scan, finds the Matching HV pairs between the two, and removes any nonmatching HV points in either one.
        moved into helpers because it is being used in multiple files
        Arguments:
            -src: hvscan with source
            -drk: hvscan without source
    '''
    src = np.array(src)
    drk = np.array(drk)
    mtchs,src_idxs,drk_idxs = np.intersect1d(src[0],drk[0],return_indices=True)
    print('\t\tMatched')
    return src[:,src_idxs],drk[:,drk_idxs]

def getPlateauMean(mscan_list, start_range = 0, end_range = 500):
    '''Helper Function used for mkGasGain and mkGasGainTable. calculates the mean value of the graph plateaus from 0(not inclusive) to 500(inclusive)'''

    mask_plateau = (mscan_list[0] > start_range) & (mscan_list[0] <= end_range)
    plateau_vals = mscan_list[1][mask_plateau]
    plateau_mean = np.mean(plateau_vals)
    plateau_std = np.std(plateau_vals)
    
    return plateau_mean,plateau_std


## To clean up  - - - - - - - - - - - - - - - - - - - - - - - - -
def mkGasGainComprehensive(src,drk,mscan_list, strip, hole,tag, start_volt=3300, end_volt=3500, start_plateau=1, end_plateau=500, x_lower_lim=None, x_upper_lim=None, y_lower_lim=None, y_upper_lim=None):    
    '''
        Create figure showing GasGain. Code largely taken from mkRawFittedPlot.

        Arguments:
            -mscan_list: list of HV, avg Current Vals, and stderr
            -strip: strip number for labeling
            -start_volt: starting voltage for  exponential fit
            -end_volt: ending voltage for exponential fit
            -start_plateau: starting point for plateau value calculation, which is used to calculate gas gain
            -end_plateau: ending point for plateau calculation
            -x_lower_lim: lower limit for graph zoom on x axis
            -x_upper_lim: upper limit for graph zoom on x axis
            -y_lower_lim: lower limit for graph zoom on y axis
            -y_upper_lim: upper limit for graph zoom on y axis

    '''
    
    #start exponential fit
    print('\n\tPerforming fit')
    mask_fit = (mscan_list[0] >= start_volt) & (mscan_list[0] <= end_volt)

    #Find mean of plateau
    plateau_mean,plateau_uncert = getPlateauMean(mscan_list, start_plateau, end_plateau)
    
    # Gas gain value
    # Gas gain = avg current / plateau mean
    gas_gain_list = mscan_list[1] / plateau_mean
    
    #get curve fit
    p0 = [0.001,0.01]
    p1,cov = curve_fit(func.mExp, mscan_list[0][mask_fit], gas_gain_list[mask_fit], p0)

    #x and y vals for exp plotting
    xs = np.linspace(0,3850,2000)
    ys = func.mExp(xs,p1[0],p1[1])

    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(src[0],src[1],src[2], marker='.',linestyle='',color='black',label='Source Current')
    plt.errorbar(drk[0],drk[1],drk[2], marker='s',markersize=3,linestyle='',color='blue',label='Dark Current')
    # plt.errorbar(mscan_list[0], mscan_list[1], marker='.',linestyle='',color='black',label='Corrected Scan')
    # plt.plot(xs,ys)

    # #show line of fit range NOTE: -1 needed to prevent error sometimes?
    # plt.axvline(start_volt,color='grey',alpha=0.3,linestyle='--')
    # plt.axvline(end_volt,color='grey',alpha=0.3,linestyle='--')
    
    #Display Plateau value used to calculate gas gain
    # plt.text(0.02, 0.98, verticalalignment='top',horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'PLATEAU VAL: {plateau_mean:.4f}', fontweight='bold', color='blue')
    
    #change scaling of plot
    # plt.xlim(x_lower_lim, x_upper_lim)
    # plt.ylim(y_lower_lim, y_upper_lim)
    # plt.yscale('symlog')

    ax = plt.gca()
    ax.set_yscale('log')
    # # ax.yaxis.set_major_formatter(mticker.FuncFormatter(sci_notation_formatter))
    # ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    lbl_fsize = 18
    mplt.add_cms_label(plt.gcf(),ax=ax,cms_fontsize=lbl_fsize,status_fontsize=lbl_fsize,context='904 Lab',context_fontsize=lbl_fsize,gap=0.1)
    # ax.text(0.03, 0.97, "CMS MUON", transform=ax.transAxes,
    #     fontsize=14, fontweight='bold', ha='left', va='top')

    # ax.text(0.03, 0.92, "Preliminary", transform=ax.transAxes,
    #     fontsize=14, style='italic', ha='left', va='top')
    
    ax.tick_params(direction='in',which='both')
    # ax.yaxis.offsetText.set_fontsize(12)
    # ax.yaxis.offsetText.set_color('gray')
    
    #Title and label graph
    # longS = strip.replace("S","Strip ")
    # plt.title(f'{longS}',fontsize=14,loc='right',x=0.95)
    # plt.title(f'{strip}, {hole.replace("_", "").replace("0","").upper()} Gas Gain',fontsize=25)

    
    plt.xlabel('High Voltage [V]',fontsize=lbl_fsize,ha='right', x=1.0)
    plt.ylabel('I [nA]',fontsize=lbl_fsize,ha='right', y=1.0)
    plt.legend()
    #Display Graph for Testing
    # plt.show()

    #Save and close graph
    # plt.tight_layout()
    plt.savefig(f'./plots/HV_Scans/raw_{strip}_{tag}.png', format='png', dpi=400)
    plt.savefig(f'./plots/HV_Scans/raw_{strip}_{tag}.pdf', format='pdf')
    plt.close()

def mkGasGain(mscan_list, strip, hole,tag, start_volt=3000, end_volt=3600, start_plateau=1, end_plateau=500, x_lower_lim=None, x_upper_lim=None, y_lower_lim=None, y_upper_lim=None):    
    '''
        Create figure showing GasGain. Code largely taken from mkRawFittedPlot.

        Arguments:
            -mscan_list: list of HV, avg Current Vals, and stderr
            -strip: strip number for labeling
            -start_volt: starting voltage for  exponential fit
            -end_volt: ending voltage for exponential fit
            -start_plateau: starting point for plateau value calculation, which is used to calculate gas gain
            -end_plateau: ending point for plateau calculation
            -x_lower_lim: lower limit for graph zoom on x axis
            -x_upper_lim: upper limit for graph zoom on x axis
            -y_lower_lim: lower limit for graph zoom on y axis
            -y_upper_lim: upper limit for graph zoom on y axis

    '''
    
    #start exponential fit
    print('\n\tPerforming fit')
    mask_fit = (mscan_list[0] >= start_volt) & (mscan_list[0] <= end_volt)

    #Find mean of plateau
    plateau_mean,plateau_std = getPlateauMean(mscan_list, start_plateau, end_plateau)
    
    # Gas gain value
    # Gas gain = avg current / plateau mean
    gas_gain_list = mscan_list[1] / plateau_mean
    
    #get curve fit
    p0 = [0.001,0.01]
    p1,cov = curve_fit(func.mExp, mscan_list[0][mask_fit], gas_gain_list[mask_fit], p0)

    #x and y vals for exp plotting
    xs = np.linspace(0,3850,2000)
    ys = func.mExp(xs,p1[0],p1[1])

    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(mscan_list[0], gas_gain_list, marker='.',linestyle='')
    plt.plot(xs,ys)

    #show line of fit range NOTE: -1 needed to prevent error sometimes?
    plt.axvline(start_volt,color='grey',alpha=0.3,linestyle='--')
    plt.axvline(end_volt,color='grey',alpha=0.3,linestyle='--')
    
    #Display Plateau value used to calculate gas gain
    # plt.text(0.02, 0.98, verticalalignment='top',horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'PLATEAU VAL: {plateau_mean:.4f}', fontweight='bold', color='blue')
    
    #change scaling of plot
    plt.xlim(x_lower_lim, x_upper_lim)
    plt.ylim(y_lower_lim, y_upper_lim)

    ax = plt.gca()
    # ax.yaxis.set_major_formatter(mticker.FuncFormatter(sci_notation_formatter))
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    mplt.add_cms_label(plt.gcf(),ax=ax,x_offset=0.081,gap=0.2,cms_fontsize=14,status_fontsize=14,context=strip.replace("S","Strip "),context_fontsize=14)
    # ax.text(0.03, 0.97, "CMS MUON", transform=ax.transAxes,
    #     fontsize=14, fontweight='bold', ha='left', va='top')

    # ax.text(0.03, 0.92, "Preliminary", transform=ax.transAxes,
    #     fontsize=14, style='italic', ha='left', va='top')
    
    ax.tick_params(direction='in')
    # ax.yaxis.offsetText.set_fontsize(12)
    # ax.yaxis.offsetText.set_color('gray')
    
    #Title and label graph
    # longS = strip.replace("S","Strip ")
    # plt.title(f'{longS}',fontsize=14,loc='right',x=0.95)
    # plt.title(f'{strip}, {hole.replace("_", "").replace("0","").upper()} Gas Gain',fontsize=25)

    
    lbl_fsize = 14
    plt.xlabel('High Voltage [V]',fontsize=lbl_fsize)
    plt.ylabel('Gas Gain',fontsize=lbl_fsize)
    #Display Graph for Testing
    # plt.show()

    #Save and close graph
    # plt.tight_layout()
    plt.savefig(f'./plots/HV_Scans/GasGain_{strip}_{tag}.png', format='png', dpi=400)
    plt.close()