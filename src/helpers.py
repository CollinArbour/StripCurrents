from os import error
from re import M
from telnetlib import XAUTH
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pyparsing import deque
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import src.DataFile as df


## Functions
def createRun(run_nm, criteria, src=True):
    '''Does Basic Data processing using the dataFile class'''
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

        Arguments:
            -src: hvscan with source
            -drk: hvscan without source
    '''
    src = np.array(src)
    drk = np.array(drk)
    mtchs,src_idxs,drk_idxs = np.intersect1d(src[0],drk[0],return_indices=True)
    print('\t\tMatched')
    return src[:,src_idxs],drk[:,drk_idxs]

def getCurveParams(mscan_list, mask_fit, p0=[0.001,0.01]):
    '''Moving curve fit calculation to this to remove space charge reliance on raw fitted'''
    p1, cov = curve_fit(mExp, mscan_list[0][mask_fit], mscan_list[1][mask_fit],p0)
    return p1

def quadSum(a,b):
    return np.sqrt(a**2 + b**2)

def Gauss(x, A, B, C):
    '''
        Returns a gaussian function given statistical info
        
        Arguments:
            @x: Independent variable values
            @A: Amplitude
            @B: Standard Deviation
            @C: Center of Distribution
    '''
    y = A*np.exp(-1/2*((x-C)/B)**2)
    return y

def mGaussianSum(x_vals,A0,B0,C0,A1,B1,C1):
    '''Create two Gaussian functions and Combine them'''
    return Gauss(x_vals,A0,B0,C0) + Gauss(x_vals,A1,B1,C1)

def mExp(x_vals,A,B):
    '''Return exponential '''
    #A*e^x
    return A*np.exp(B*x_vals)

def m2DGaussian(x_vals,y_vals,A0,B0,C0,A1,B1,C1):
    '''Creates two 2d gaussians and combines them'''
    return A0*Gauss(x_vals,1,B0,C0)*Gauss(y_vals,1,B0,C0) + A1*Gauss(x_vals,1,B1,C1)*Gauss(y_vals,1,B1,C1)

def FWHM(xs,ys):
    '''Returns the x-values of which the full width half max lies'''
    hm = np.max(ys)/2
    x1 = 0
    x1d = 100
    x2 = 0
    x2d = 100

    for i,y in enumerate(ys):
        d = np.abs(y-hm)

        if xs[i] < 76:
            if d < x1d:
                x1 = xs[i]
                x1d = d
        else:
            if d < x2d:
                x2 = xs[i]
                x2d = d
    
    return x1,x2

def intFWHM(ps,fwhm,pts=1000):
    '''
        Function should probably be phased out in favor of 'intRadius(...)'
    '''
    a0 = ps[0] / (np.sqrt(2*np.pi) * ps[1]) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = ps[3] / (np.sqrt(2*np.pi) * ps[4]) # Amplitude of second 2D Gaussian (nA/mm^2)
    
    r = (fwhm[1]-fwhm[0])/2
    
    G0,G0_err = integrate.quad(lambda x: Gauss(x,1,ps[1],0),-r,r)
    G1,G1_err = integrate.quad(lambda x: Gauss(x,1,ps[4],0),-r,r)
    
    G2D = a0*G0**2 + a1*G1**2
    
    return G2D

def intRadius(ps,r,pts=1000):
    '''
    Integrates total current within area of circle r

    Arguments:
        @ps : parameters of sum of two 2D Gaussians [a0,sig0,0,a1,sig1,0]
        @r  : Radius of circle to perform  integration over
        @pts : Number of points to integrate over
    '''
    
    G0,G0_err = integrate.quad(lambda x: Gauss(x,1,ps[1],0),-r,r)
    G1,G1_err = integrate.quad(lambda x: Gauss(x,1,ps[4],0),-r,r)
    
    G2D = ps[0]*G0**2 + ps[3]*G1**2
    
    return G2D

def intRadiusCylindrical(ps,r,pts=1000):
    '''
    Integrate in cylindrical coordinates
    '''
    G0,G0_err = integrate.quad(lambda x: x*Gauss(x,1,ps[1],0),0,r)
    G1,G1_err = integrate.quad(lambda x: x*Gauss(x,1,ps[4],0),0,r)
    
    G2D = ps[0] *(2*np.pi)*G0 + ps[3] *(2*np.pi)*G1

    return G2D

def wirePlacement(r):
    wireGap = 3.13 #mm
    wireD = 0.050 #mm
    wireUnit = wireGap + wireD
    wire_xs = np.arange(-r,r,wireUnit)

    return wire_xs

def wireLength(r,offset=None):
    xs = wirePlacement(r)
    ltot = 0
    for x in xs:
        ltot += 2*np.sqrt(r**2-x**2)

    return ltot

def accumCharge(Itot,r):
    '''
    Returns accumulated  charge in (mC/cm) / day

    Arguments
        @Itot : Total current in beamspot of radius r in nA
        @r :  radius of beam spot
    '''
    lwire = wireLength(r)

    Imc = 2*Itot/1000000 # nA -> mA
    lcm = lwire/10 #mm -> cm

    return Imc/lcm * 60 * 60 * 24
    
def mkScans(strips,ps,i,save=False):
    '''
    arg:strips - [xpos,values,stderrs] **values is nA/mm (across with width of the stri)
    arg:ps - curve fit
    '''
    # Producing line shape of fit
    xfit = np.linspace(0,200,1000)
    yfit = mGaussianSum(xfit,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5])
    fwhm = FWHM(xfit,yfit)
    fwhmCharge = intFWHM(ps,fwhm)
    
    plt.errorbar(strips[0],strips[1],yerr=strips[2],linestyle='--')
    plt.plot(xfit,yfit) # ,label=nms[i]) <--- NEED TO ADAPT FUNCTION TO TAKE IN LABELS
#     for x in fwhm:
#         plt.axvline(x,linestyle=':',color='grey',alpha=0.5)
    
    plt.xlabel('Distance (mm)')
    plt.ylabel('Linear Current Density (nA/mm)')
    plt.title(f'MiniCSC 90Sr L1 H2 HV3600, Strip Scan')
    # plt.title(f'MiniCSC 90Sr Src{i+1} L1 H2 HV3600, Strip Scan')
    plt.legend()
    
    plt.text(110,23-i*3.5,f'Src{i+1} FWHM Itot: {fwhmCharge:.2f} nA',bbox=dict(facecolor='grey',alpha=0.75))
    
    if save:
        plt.savefig(f'./plots/SrSrcs/avgI_src{i}.png',format='png',dpi=400)
        plt.close()
    else:
        plt.show()

def mkHeatMap_GaussSum(r,ps,pts=1000,mlabel='',save=False):
    '''
    Makes heat map showing 2D distribution of current from source

    Arguments:
        @r : Distance from center point of source that map will cover
        @ps : Parameters of fit to sum of two Gaussians
        @pts : Number of points between -r and r to use when plotting
        @mlabel : label
    '''
    sig0 = ps[1]
    sig1 = ps[4]
    a0 = ps[0] / (np.sqrt(2*np.pi) * sig0) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = ps[3] / (np.sqrt(2*np.pi) * sig1) # Amplitude of second 2D Gaussian (nA/mm^2)

    # Determine FWHM
    xfit = np.linspace(0,200,1000)
    yfit = mGaussianSum(xfit,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5])
    fwhm = FWHM(xfit,yfit)
    print(fwhm)
    fwhm_d = fwhm[1]-fwhm[0]
    fwhm_r = fwhm_d/2

    #creates grid
    hmpts = np.linspace(-r,r,pts)
    hmxpts,hmypts = np.meshgrid(hmpts,hmpts)
    
    G0x = Gauss(hmxpts,1,sig0,0)
    G0y = Gauss(hmypts,1,sig0,0)
    
    G1x = Gauss(hmxpts,1,sig1,0)
    G1y = Gauss(hmypts,1,sig1,0)
    
    G = a0*G0x*G0y + a1*G1x*G1y
    
    hm = np.max(G) / 2

    sig0_lvl = 5        
    sig0_lvl_xval = sig0 * sig0_lvl     #5sigma
    sig0_lvl_yval = mGaussianSum(sig0_lvl_xval,a0,sig0,0,a1,sig1,0)     #Create and sum 2 Gaussians 

    contour = plt.contour(hmxpts, hmypts, G, levels=[sig0_lvl_yval,hm], colors='white', linestyles='dashed', linewidths=2)
    plt.clabel(contour, inline=True, fontsize=8, fmt={sig0_lvl_yval: f'{sig0_lvl}σ_0',hm: 'FWHM'})

    sig0_lvl_totChrg = intRadius([a0,sig0,0,a1,sig1,0],sig0_lvl_xval)   #5sigma
    sig0_lvl3_totChrg = intRadius([a0,sig0,0,a1,sig1,0],sig0*3)         #3sigma
    sig0_lvl10_totChrg = intRadius([a0,sig0,0,a1,sig1,0],sig0*10)       #10sigma
    
    if np.max(G) > 1.5:
        print('ERROR: Max is greater than upper limit in Heat Map')

    plt.imshow(G,origin='lower', cmap='viridis',extent=[-r,r,-r,r],vmin=0,vmax=1.5)
    
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')
    plt.title(f'Sum of 2D-Gaussians Fit {mlabel}')

    xtxtpos = -41
    ytxtpos = -41
    plt.text(xtxtpos,ytxtpos,f'σ0: {sig0:0.2f}mm  σ1: {sig1:.2f}mm  FWHM: \nItot in 3*σ0= {sig0_lvl3_totChrg:.2f} nA\nItot in 5*σ0= {sig0_lvl_totChrg:.2f} nA\nItot in 10*σ0= {sig0_lvl10_totChrg:.2f} nA',bbox=dict(facecolor='grey',alpha=1))
    #plt.text(xtxtpos,ytxtpos,f'r=5σ_0={sig:.2f} mm \nItot: {chrg_sig:.2f} nA',bbox=dict(facecolor='grey',alpha=0.75))
    
    cbar = plt.colorbar()  # Create a colorbar
    cbar.set_label('Current Density (nA/mm^2)', rotation=270, labelpad=15)
    
    if mlabel and save:
        mlabel = '_'.join(mlabel.split())
        plt.savefig(f'./plots/SrSrcs/{mlabel}_HeatMap.png',format='png',dpi=400)
        plt.close()
    else:
        plt.show()

def mkGasGain(mscan_list, strip, hole, start_volt=3000, end_volt=3600, start_plateau=0, end_plateau=500, x_lower_lim=None, x_upper_lim=None, y_lower_lim=None, y_upper_lim=None, plot=True):    
    '''
        Create figure showing GasGain. Code largely taken from mkRawFittedPlot.
        returns plateau_mean and gas gain val for making table(for now)

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
            -plot: whether to plot, might not even implement this and instead just move calculations to another method so this doesnt have to be called.

    '''
    
    #start exponential fit
    print('\n\tPerforming fit')
    mask_fit = (mscan_list[0] >= start_volt) & (mscan_list[0] <= end_volt)

    #Find index of plateau starting and ending points
    mask_plateau = (mscan_list[0] > start_plateau) & (mscan_list[0] <= end_plateau)

    #Find mean of plateau
    plateau_vals = mscan_list[1][mask_plateau]
    plateau_mean = np.mean(plateau_vals) 
    
    # Gas gain value
    # Gas gain = avg current / plateau mean
    mscan_list[1] /= plateau_mean
    
    #get curve fit
    p0 = [0.001,0.01]
    print("TEST HV RANGE: ", mscan_list[0][mask_fit])
    print("TEST GASGAIN RANGE: ", mscan_list[1][mask_fit])

    p1,cov = curve_fit(mExp, mscan_list[0][mask_fit], mscan_list[1][mask_fit], p0)

    #x and y vals for exp plotting
    xs = np.linspace(0,3850,2000)
    ys = mExp(xs,p1[0],p1[1])

    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(mscan_list[0], mscan_list[1], marker='.',linestyle='')
    plt.plot(xs,ys)

    #show line of fit range NOTE: -1 needed to prevent error sometimes?
    plt.axvline(start_volt,color='grey',alpha=0.3,linestyle='--')
    plt.axvline(end_volt,color='grey',alpha=0.3,linestyle='--')
    
    #Display Plateau value used to calculate gas gain
    plt.text(0.02, 0.98, verticalalignment='top',horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'PLATEAU VAL: {plateau_mean:.4f}', fontweight='bold', color='blue')
    
    #change scaling of plot
    plt.xlim(x_lower_lim, x_upper_lim)
    plt.ylim(y_lower_lim, y_upper_lim)
    
    #Title and label graph
    plt.title(f'{strip}, {hole.replace("_", "").replace("0","").upper()} Gas Gain')
    plt.xlabel('HV (V)')

    plt.ylabel('Gas Gain')
    #Display Graph for Testing
    #plt.show()

    #Save and close graph
    #plt.savefig(f'./plots/HV_Scans/GasGain/Graphs/{strip}.png', format='png', dpi=400)
    plt.close()

    combined_points = zip(mscan_list[0], mscan_list[1])
   
    gas_gain_val = list(point[1] for point in combined_points if point[0] == 3600)
    return plateau_mean, gas_gain_val
    
def mkRawFittedPlot(mscan_list, strip, start_volt=3000, end_volt=3550):
    '''
        This Function serves 2 Purposes:
            -Plots raw fitted graph
            -Returns p1 parameter list
                - p1 is necessary to run mkSpaceChargePlot

        Arguments: 
            mscan_list: list that includes HV points, Corrected Avg current points, And stderr. this is just mhvscan in mkHVScanPlot.
            strip: which strip is being irradiated, for graph labeling
            plot: determines whether to generate the plot, or just return p1.

    '''
    print('\n\tPerforming fit')
    mask_fit = (mscan_list[0] >= start_volt) & (mscan_list[0] <= end_volt)
    p1 = getCurveParams(mscan_list, mask_fit)
    
    xs = np.linspace(0,3850,2000)
    ys = mExp(xs,p1[0],p1[1])
    # Create raw figure
    print('\tCreating raw figure')
    plt.errorbar(mscan_list[0],mscan_list[1],yerr=mscan_list[2],marker='.',linestyle='')
    plt.plot(xs,ys)
    plt.axvline(start_volt, color='grey',alpha=0.3,linestyle='--')
    plt.axvline(end_volt, color='grey',alpha=0.3,linestyle='--')
    plt.title(f'{strip} Strip Current over HV Scan')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.show()
    #plt.savefig(f'./plots/HV_Scans/{strip}_RawFitted.png',format='png',dpi=400)
    plt.close()

def mkPlateauPlot(mscan_list, strip, src, hole, **kwargs):       
    '''
        This function produces a graph of the plateau, including error bars and a linear fit. 
        It also calculates the average value of the plateau and labels it.

        Arguments:
            mscan_list: list that includes HV points, Corrected Avg current points, And stderr. this is just mhvscan in mkHVScanPlot.
            strip: which strip is being irradiated, for graph labeling
            src: source used during measurements, for graph labeling
            hole: hole used during measurements, for graph labeling

        Optional Arguments:
            start_point: starting hv value, defaulted to zero
            end_point: ending hv value, defaulted to 500.
            lim args: to manually change the zoom on the graph
            
            uncorrected_curr: takes in the src_hvscan in mkHVScanPlot to graph the uncorrected source current
            uncorrected_dark_curr: takes in the drk_hvscan in mkHVScanPlot to graph the uncorrected source current
                -NOTE: for this and the parameter above it:
                    -index 0 holds the HV (x-axis) values
                    -index 1 holds the Current (y-axis) values

        
    '''
    defaults = {
        #Dict holds more parameters. accessed by  value = defaults['key']

        #optional parameters
        'end_point': 500,           #(V), ending point in data
        'start_point': 5,           #(V), starting point in data run.
        #change scale of graph image
        'x_low_lim': None,          
        'x_upper_lim': None,
        'y_lower_lim': None,
        'y_upper_lim': None,
        #to graph pre-corrected src and dark values
        'uncorrected_curr': None,
        'uncorrected_dark_curr': None
    }
    #set parameters in
    defaults.update(kwargs)
    
    #notify user of plateau figure cpreation
    print("Creating Plateau Figure")
    
    mask_fit = (mscan_list[0] >= defaults["start_point"]) & (mscan_list[0] <= defaults['end_point'])

    #grab list of plateau avg curr vals 
    plateau_vals = mscan_list[1][mask_fit]
    plateau_mean = np.mean(plateau_vals)

    #grab x and y values of points
    x_vals = mscan_list[0][mask_fit]
    y_vals = mscan_list[1][mask_fit]

    #plot value points with errorbars
    plt.errorbar(x_vals, y_vals,yerr=mscan_list[2][mask_fit], marker='.', linestyle='', label='Corrected Avg Current')
    
    
    if (defaults['uncorrected_curr'] is not None) and (defaults['uncorrected_dark_curr'] is not None):
        
        #grab x and y values of uncorrected source points
        x_uncorrected = defaults['uncorrected_curr'][0][mask_fit]      #first bracket accesses defaults dict value, 2nd bracket accesses HV values(for x-axis), 3rd bracket accesses what range of values are needed
        y_uncorrected = defaults['uncorrected_curr'][1][mask_fit]      #first bracket accesses defaults dict value, 2nd bracket accesses Avg Curr values(for y-axis), 3rd bracket accesses what range of values are needed

        #grab x and y values of uncorrected dark points
        x_uncorrected_dark = defaults['uncorrected_dark_curr'][0][mask_fit]        #Check comments above^^
        y_uncorrected_dark = defaults['uncorrected_dark_curr'][1][mask_fit]        #Check comments above^^

        #grab the standard error of the uncorrected runs for the errorbars
        src_bar_yerr = abs(defaults['uncorrected_curr'][2][mask_fit])
        drk_bar_yerr = abs(defaults['uncorrected_dark_curr'][2][mask_fit])

        #plot uncorrected src and drk dat points
        plt.errorbar(x_uncorrected, y_uncorrected, yerr = src_bar_yerr, marker='.', linestyle='', color='green', label='Uncorrected Src Avg Current')
        plt.errorbar(x_uncorrected_dark, y_uncorrected_dark, yerr = drk_bar_yerr, marker='.',linestyle='', color='orange', label='Uncorrected Dark Avg Current')

    #add text showing the mean of the plateau
    plt.text(0.02, 0.98, verticalalignment='top',horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'{defaults["start_point"]:.0f}V-{defaults["end_point"]:.0f}V: {plateau_mean:.4f} nA', fontweight='bold', color='blue')
    
    #calculate linear fit
    slope, intercept = np.polyfit(x_vals, y_vals, 1)
    fit_line = slope * x_vals + intercept

    #plot linear fit
    plt.scatter(x_vals, y_vals)
    plt.plot(x_vals, fit_line, color='red')

    #manually change zoom of graph, automatically done if passed all "None" values
    plt.xlim(defaults['x_low_lim'], defaults['x_upper_lim'])
    plt.ylim(defaults['y_lower_lim'], defaults['y_upper_lim'])

    #add text of Linear Fit Function
    plt.text(0.02, 0.9, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'y = {slope:.2e} * x + {intercept:.2e}', fontweight='bold', color='black')

    #label and title plot, save and close
    plt.title(f'MiniCSC4: HV Scan, L1, 90{src}-Src1, {hole.replace("_", "").replace("0","").upper()}, {strip}, Space Charge Fit')
    plt.xlabel('HV (V)')
    plt.ylabel('ln(Avg. I)')
    plt.legend()
    plt.show()
    #plt.savefig(f'./plots/HV_Scans/{strip}_src{src}_plateau.png', format='png', dpi=400)
    plt.close()

def mkSpaceChargePlot(mscan_list, strip, start_volt=3000, end_volt=3600):
    '''
        Creates Space Charge graph ****needs more explanation

        Arguments:
            mscan_list: list that includes HV points, Corrected Avg current points, And stderr. this is just mhvscan in mkHVScanPlot.
            strip: which strip is being irradiated, for graph labeling
    '''

    fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 3]}, sharex=True)
    fig.subplots_adjust(hspace=0)

    mask_lower = (mscan_list[0] >= start_volt)
    mask_upper = (mscan_list[0] <= end_volt)

    p1 = getCurveParams(mscan_list, (mask_lower * mask_upper),)

    xs = np.linspace(3000, 3850, 1500)
    ys = mExp(xs, p1[0], p1[1])
     
    #TODO: Figure out how to access the beginning of the mask but not cutoff(chatgpt)
    print('\tCreating Space Charge evaluation')
    ax0.errorbar(mscan_list[0][mask_lower],mscan_list[1][mask_lower],yerr=mscan_list[2][mask_lower],marker='.',linestyle='')
    ax0.plot(xs,ys)
    ax0.set_yscale('symlog')
    ax0.set_title(f'{strip} Strip Current over HV Scan')
    ax0.set_ylabel('Avg. I (nA)')

    yexpect = mExp(mscan_list[0][mask_lower],p1[0],p1[1])
    diff = yexpect - mscan_list[1][mask_lower]
    rel_diff = diff / mscan_list[1][mask_lower]

    ax1.plot(mscan_list[0][mask_lower],rel_diff)
    ax1.set_xlabel('HV (V)')
    ax1.set_ylabel('(Exp-Data)/Data')

    ax0.axvline(start_volt,color='grey',alpha=0.3,linestyle='--')
    ax0.axvline(end_volt,color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(start_volt,color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(end_volt,color='grey',alpha=0.3,linestyle='--')
    ax1.axhline(0,color='black',alpha=0.5,linestyle=':')

    plt.show()
    #plt.savefig(f'./plots/HV_Scans/{strip}_SpaceCharge.png',format='png',dpi=400)
    plt.close()

def mkGasGainTable(table_plateaus, table_gasGain, table_hole, table_strip):
    '''
        Creates a table of values to show gas gain across different runs.

        Arguments:
            -table_plateaus: numpy array of plateau mean values (returned from mkGasGain and put in as parameter)
            -table_gasGain: numpy array of gasGain value at 3600v for each run (also returned from mkGasGain and put in as parameter)
            -table_holes: array of hole values for labeling holes
            -table_strip: array of strip values for labeling strips
    '''

    #Round the plateau mean values and the gas gain values
    rounded_plateaus = np.round(table_plateaus, 4)
    rounded_gasGain = np.round(table_gasGain, 1)

    #removes brackets from gasGain values, for cleaner table visual
    formatted_gasGain = []
    for val in rounded_gasGain:
        formatted_gasGain.append(val[0])

    #adds commas for gas gain values, for readability
    formatted_gasGain = np.vectorize(lambda x: f"{x:,}")(formatted_gasGain)

    #creates the row labels, showing both hole and strip
    formatted_labels = []
    for hole, strip in zip(table_hole, table_strip):
        formatted_labels.append(hole.replace("_", "").replace("0","").upper() + ', ' + strip)
        
    #Grab inputs for table creation
    data_rows = list(zip(rounded_plateaus, formatted_gasGain))  #zip to make row of data cells, then convert to list
    row_labels = formatted_labels                               #labels for each row
    column_labels = ['Plateau Avg', 'GasGain Values']                  #labels for each column

    #formatting. sets size, and removes graph things to make only table
    fig, ax = plt.subplots(figsize=(5, 2))
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.axis('off')

    #make table
    table = ax.table(cellText=data_rows, colLabels=column_labels, rowLabels=row_labels, cellLoc='center', loc='center')

    #label and scale table
    plt.title('GasGain Table (Values at 3600V)')
    table.scale(1, 1.5)

    #Formatting, removes weird box outline
    for key, cell in table.get_celld().items():
        if key[0] == -1 or key[1] == -1:  # Header or index row/column
            cell.set_linewidth(1)
        else:
            cell.set_linewidth(0.5)

    plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  
    

    #show table for testing
    #plt.show()
    
    #Save and close table
    plt.savefig(f'./plots/HV_Scans/GasGain/Tables/GasGain_refMeasures.png',format='png',dpi=400)
    plt.close()
