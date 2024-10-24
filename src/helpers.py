import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as integrate


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
    

def mkScans(strips,ps,i,save=False,markers=True):
    '''
    @arg:strips - [xpos,values,stderrs]
    @arg:ps - []
    '''
    
    # a0 = ps[0] / (np.sqrt(2*np.pi) * ps[1]) # Amplitude of first 2D Gaussian (nA/mm^2)
    # a1 = ps[3] / (np.sqrt(2*np.pi) * ps[4]) # # Amplitude of second 2D Gaussian (nA/mm^2)
    
    # Producing line shape of fit
    xfit = np.linspace(0,200,1000)
    yfit = mGaussianSum(xfit,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5])
    fwhm = FWHM(xfit,yfit)
    fwhmCharge = intFWHM(ps,fwhm)
    
    # yfit2 = mGaussianSum(xfit,a0,ps[1],ps[2],a1,ps[4],ps[5])
    # creating the graph
    # zorder is used si that the error bars show on top of the function
    #
    fig, ax = plt.subplots()

    if markers:
        ax.errorbar(strips[0],strips[1],yerr=strips[2],linestyle=' ', marker='.', zorder=3)
        ax.text(0.03,0.93,f'Primary',bbox=dict(facecolor='white',alpha=0.5),transform=ax.transAxes)

    else:    
        ax.text(0.03,0.93,f'σ(x, y=0)',bbox=dict(facecolor='white',alpha=0.5),transform=ax.transAxes)
    
    # ax.text(0.6,0.725,f'Src{i+1} FWHM Itot: {fwhmCharge:.2f} nA',bbox=dict(facecolor='grey',alpha=0.75),transform=ax.transAxes)
    
    ax.plot(xfit,yfit, color='black', alpha=0.75, zorder=2) 
    ax.axis([0, max(xfit), 0, max(yfit) * 1.1])
    
    # ,label=nms[i]) <--- NEED TO ADAPT FUNCTION TO TAKE IN LABELS
#     for x in fwhm:
#         plt.axvline(x,linestyle=':',color='grey',alpha=0.5)
    
    ax.set_xlabel('Distance (mm)')
    ax.set_ylabel('Linear Current Density (nA/mm)')
    ax.set_title(f'MiniCSC4: Strip Scan L1 90Sr-Src{i+1} H2 HV3600')
    '''
    title format:

    MiniCSC4 {Measurement type} L1 H2 90Sr-Src{i+1} HV3600
            strip scan vs hvSCAN 
    '''
    
    # plt.title(f'MiniCSC 90Sr Src{i+1} L1 H2 HV3600, Strip Scan')
    

    # Putting denotations on the graph for FWHM, 1Sig and 2Sig
    # ps[2] is the mean
    # ps[1] is sigma
    mu = ps[2]
    sigma = ps[1]
    
    ax.fill_between(xfit, yfit, where=((xfit >= fwhm[0]) & (xfit <= fwhm[1])), color='red', alpha=0.4, label='FWHM range')
    ax.fill_between(xfit, yfit, where=((xfit >= mu - 2 * sigma) & (xfit <= mu + 2 * sigma)), color='red', alpha=0.3, label='2σ range')
    ax.fill_between(xfit, yfit, where=((xfit >= mu - 3 * sigma) & (xfit <= mu + 3 * sigma)), color='red', alpha=0.2, label='3σ range')
    
    ax.legend()
    
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
