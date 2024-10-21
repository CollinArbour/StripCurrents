import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as integrate
from scipy.optimize import curve_fit


## Functions
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
    @arg:strips - [xpos,values,stderrs]
    @arg:ps - []
    '''
    # Producing line shape of fit
    xfit = np.linspace(0,200,1000)
    yfit = mGaussianSum(xfit,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5])
    fwhm = FWHM(xfit,yfit)
    fwhmCharge = intFWHM(ps,fwhm)
    
    plt.errorbar(strips[0],strips[1],yerr=strips[2],linestyle='',marker='.')
    plt.plot(xfit,yfit) # ,label=nms[i]) <--- NEED TO ADAPT FUNCTION TO TAKE IN LABELS
    
    plt.xlabel('Distance (mm)')
    plt.ylabel('Linear Current Density (nA/mm)')
    #plt.ylabel('Avg. Current (nA)')
    plt.title(f'MiniCSC4: Strip Scan L1 90Sr-Src{i+1} H2 HV3600')
    plt.legend()
    
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

    sig0_lvl = 3
    sig0_lvl_xval = sig0 * sig0_lvl
    sig0_lvl_yval = mGaussianSum(sig0_lvl_xval,a0,sig0,0,a1,sig1,0)

    sig0_lvl01 = 2
    sig0_lvl_xval01 = sig0 * sig0_lvl01
    sig0_lvl_yval01 = mGaussianSum(sig0_lvl_xval01,a0,sig0,0,a1,sig1,0)

    contour = plt.contour(hmxpts, hmypts, G, levels=[sig0_lvl_yval,sig0_lvl_yval01,hm], colors='white', linestyles='dashed', linewidths=2)
    plt.clabel(contour, inline=True, fontsize=8, fmt={sig0_lvl_yval: f'{sig0_lvl}σ_0',sig0_lvl_yval01: f'{sig0_lvl01}σ_0',hm: 'FWHM'})

    sig0_lvl_totChrg = intRadiusCylindrical([a0,sig0,0,a1,sig1,0],sig0_lvl_xval)
    sig0_lvl01_totChrg = intRadiusCylindrical([a0,sig0,0,a1,sig1,0],sig0_lvl_xval01)
    fwhm_totChrg = intRadiusCylindrical([a0,sig0,0,a1,sig1,0],fwhm_r)
    
    if np.max(G) > 1.5:
        print('ERROR: Max is greater than upper limit in Heat Map')

    #mnorm = colors.SymLogNorm(linthresh=0.01,vmin=0,vmax=1.5)
    #plt.imshow(G,origin='lower', cmap='viridis',norm=mnorm,extent=[-r,r,-r,r]) #,vmin=0,vmax=1.5)
    plt.imshow(G,origin='lower', cmap='viridis',extent=[-r,r,-r,r]) #,vmin=0,vmax=1.5)
    
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')
    plt.title(f'Sum of 2D-Gaussians Fit {mlabel}')

    xtxtpos = -35
    ytxtpos = 24
    plt.text(xtxtpos,ytxtpos,f'σ0: {sig0:0.2f}mm  σ1: {sig1:.2f}mm FWHM: {fwhm_d:0.2f}mm\nItot in FWHM= {fwhm_totChrg:.2f} nA\nItot in {sig0_lvl01}*σ0= {sig0_lvl01_totChrg:.2f} nA\nItot in {sig0_lvl}*σ0= {sig0_lvl_totChrg:.2f} nA',bbox=dict(facecolor='grey',alpha=1))
    #plt.text(xtxtpos,ytxtpos,f'r=5σ_0={sig:.2f} mm \nItot: {chrg_sig:.2f} nA',bbox=dict(facecolor='grey',alpha=0.75))
    
    cbar = plt.colorbar()  # Create a colorbar
    cbar.set_label('Current Density (nA/mm^2)', rotation=270, labelpad=15)
    
    if mlabel and save:
        mlabel = '_'.join(mlabel.split())
        plt.savefig(f'./plots/SrSrcs/{mlabel}_HeatMap.png',format='png',dpi=400)
        plt.close()
    else:
        plt.show()

def mkRawFittedPlot(mscan_list, strip, plot=True, start_volt=3000, end_volt=3550):
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
    strt_fit = np.where(mscan_list[0]==start_volt)[0][0]
    stop_fit = np.where(mscan_list[0]==end_volt)[0][0]+1
    p0 = [0.001,0.01]
    p1,cov = curve_fit(mExp, mscan_list[0][strt_fit:stop_fit], mscan_list[1][strt_fit:stop_fit],p0)
    if plot:
        xs = np.linspace(0,3850,2000)
        ys = mExp(xs,p1[0],p1[1])

        # Create raw figure
        print('\tCreating raw figure')
        plt.errorbar(mscan_list[0],mscan_list[1],yerr=mscan_list[2],marker='.',linestyle='')
        plt.plot(xs,ys)

        plt.axvline(mscan_list[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
        plt.axvline(mscan_list[0][stop_fit],color='grey',alpha=0.3,linestyle='--')

        plt.title(f'{strip} Strip Current over HV Scan')
        plt.xlabel('HV (V)')
        plt.ylabel('Avg. I (nA)')

        plt.savefig(f'./plots/HV_Scans/{strip}_RawFitted.png',format='png',dpi=400)
        plt.close()
    return p1

def mkPlateauPlot(mscan_list, strip, src, hole, **kwargs):       
    '''
        This function produces a graph of the plateau, including error bars and a linear fit. 
        It also calculates the average value of the plateau and labels it.

        Arguments:
            mscan_list: list that includes HV points, Corrected Avg current points, And stderr. this is just mhvscan in mkHVScanPlot.
            strip: which strip is being irradiated, for graph labeling
            src: source used during measurements, for graph labeling

        Optional Arguments:
            num_calcs: number of times to perform a mean current calculation. currently not implemented.
            start_point: starting hv value, defaulted to zero
            exclude_start: whether or not to exclude the first point(mainly to exclude the 0V run). should probably be changed to false if modifying start_point.
            end_point: ending hv value, defaulted to 500.
            lim args: to manually change the zoom on the graph
            
            uncorrected_curr: takes in the src_hvscan in mkHVScanPlot to graph the uncorrected source current
            uncorrected_dark_curr: takes in the drk_hvscan in mkHVScanPlot to graph the uncorrected source current
                -NOTE: for this and the parameter above it:
                    -index 0 holds the HV (x-axis) values
                    -index 1 holds the Current (y-axis) values

        TODO:
            -implement the num_calcs feature
                -iterative system to calculate plateau mean across different ranges, and label the value
                    -more parameters needed for range points
    '''
    defaults = {
        '''
            This Dictionary holds more parameters(all optional) using the kwargs feature of python.
            Access values like any other Dictionary
                -EX: var = defaults['num_calcs']
            When calling the function, these parameters can be changed the same way as any others. 
        '''
        #optional parameters
        'num_calcs': 1,             #number of mean calculations to do.yet to be implemented, so keep at 1.
        'start_point': 0,           #(V), starting point in data run.
        'exclude_start' : True,     #exludes starting point of data(intended to remove zero)
        'end_point': 500,           #(V), ending point in data
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
    
    #notify user of plateau figure creation
    print("Creating Plateau Figure")
    
    #Starting and ending fit points
    strt_fit = np.where(mscan_list[0] == defaults['start_point'])[0][0]   
    stop_fit = np.where(mscan_list[0] == defaults['end_point'])[0][0] + 1  #+1 to include point

    #add one to exclude starting point
    strt_fit += 1 if defaults['exclude_start'] == True else None

    #grab list of plateau avg curr vals 
    plateau_vals = mscan_list[1][strt_fit:stop_fit]
    plateau_mean = np.mean(plateau_vals)

    #grab x and y values of points
    x_vals = mscan_list[0][strt_fit:stop_fit]
    y_vals = mscan_list[1][strt_fit:stop_fit]

    #plot value points with errorbars
    plt.errorbar(x_vals, y_vals, yerr=mscan_list[2][strt_fit:stop_fit], marker='.', linestyle='', label='Corrected Avg Current')

    
    if (defaults['uncorrected_curr'] is not None) and (defaults['uncorrected_dark_curr'] is not None):
        
        #grab x and y values of uncorrected source points
        x_uncorrected = defaults['uncorrected_curr'][0][strt_fit:stop_fit]      #first bracket accesses defaults dict value, 2nd bracket accesses HV values(for x-axis), 3rd bracket accesses what range of values are needed
        y_uncorrected = defaults['uncorrected_curr'][1][strt_fit:stop_fit]      #first bracket accesses defaults dict value, 2nd bracket accesses Avg Curr values(for y-axis), 3rd bracket accesses what range of values are needed

        #grab x and y values of uncorrected dark points
        x_uncorrected_dark = defaults['uncorrected_dark_curr'][0][strt_fit:stop_fit]        #Check comments above^^
        y_uncorrected_dark = defaults['uncorrected_dark_curr'][1][strt_fit:stop_fit]        #Check comments above^^

        #grab the standard error of the uncorrected runs for the errorbars
        src_bar_yerr = abs(defaults['uncorrected_curr'][2][strt_fit:stop_fit])
        drk_bar_yerr = abs(defaults['uncorrected_dark_curr'][2][strt_fit:stop_fit])

        #plot uncorrected src and drk dat points
        plt.errorbar(x_uncorrected, y_uncorrected, yerr = src_bar_yerr, marker='.', linestyle='', color='green', label='Uncorrected Src Avg Current')
        plt.errorbar(x_uncorrected_dark, y_uncorrected_dark, yerr = drk_bar_yerr, marker='.',linestyle='', color='orange', label='Uncorrected Dark Avg Current')

    #add text showing the mean of the plateau
    plt.text(0.02, 0.98, verticalalignment='top',horizontalalignment='left', bbox=dict(facecolor='lightblue', alpha=0.5), transform=plt.gca().transAxes, s=f'{mscan_list[0][strt_fit]:.0f}V-{mscan_list[0][stop_fit]:.0f}V: {plateau_mean:.4f} nA', fontweight='bold', color='blue')
    
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
    plt.title(f'MiniCSC4: HV Scan, L1, 90{src}-Src1, {hole.replace("_", "").replace("0","").upper()}, {strip}')
    plt.xlabel('HV (V)')
    plt.ylabel('Avg. I (nA)')
    plt.legend()
    plt.savefig(f'./plots/HV_Scans/{strip}_src{src}_plateau.png', format='png', dpi=400)
    plt.close()

def mkSpaceChargePlot(mscan_list, strip, p1, start_volt=3000, end_volt=3550):
    '''
        Creates Space Charge graph ****needs more explanation

        Arguments:
            mscan_list: list that includes HV points, Corrected Avg current points, And stderr. this is just mhvscan in mkHVScanPlot.
            strip: which strip is being irradiated, for graph labeling
            p1: set of parameter values
                -THIS IS RETRIEVED FROM mkRawFittedPlot(). that code will return p1 AND print a graph. if no graph is needed, pass it plot=False in the function call
                -Limitation: Must use the same start and end voltages as the mkRawFittedPlot. 
    '''
    fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 3]}, sharex=True)
    fig.subplots_adjust(hspace=0)

    xs = np.linspace(3000, 3850, 1500)
    ys = mExp(xs, p1[0], p1[1])

    strt_fit = np.where(mscan_list[0] == start_volt)[0][0]
    stop_fit = np.where(mscan_list[0] == end_volt)[0][0] + 1

    print('\tCreating Space Charge evaluation')
    ax0.errorbar(mscan_list[0][strt_fit:],mscan_list[1][strt_fit:],yerr=mscan_list[2][strt_fit:],marker='.',linestyle='')
    ax0.plot(xs,ys)
    ax0.set_yscale('symlog')
    ax0.set_title(f'{strip} Strip Current over HV Scan')
    ax0.set_ylabel('Avg. I (nA)')

    yexpect = mExp(mscan_list[0][strt_fit:],p1[0],p1[1])
    diff = yexpect - mscan_list[1][strt_fit:]
    rel_diff = diff / mscan_list[1][strt_fit:]

    ax1.plot(mscan_list[0][strt_fit:],rel_diff)
    ax1.set_xlabel('HV (V)')
    ax1.set_ylabel('(Exp-Data)/Data')

    ax0.axvline(mscan_list[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    ax0.axvline(mscan_list[0][stop_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(mscan_list[0][strt_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axvline(mscan_list[0][stop_fit],color='grey',alpha=0.3,linestyle='--')
    ax1.axhline(0,color='black',alpha=0.5,linestyle=':')

    plt.savefig(f'./plots/HV_Scans/{strip}_SpaceCharge.png',format='png',dpi=400)
    plt.close()

