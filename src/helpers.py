import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


## Functions
def Gauss(x, A, B, C):
    '''
        A - Amplitude
        B - Standard Deviation
        C - Center of Distribution
    '''
    y = A*np.exp(-1/2*((x-C)/B)**2)
    return y

def mGaussianSum(x_vals,A0,B0,C0,A1,B1,C1):
    return Gauss(x_vals,A0,B0,C0) + Gauss(x_vals,A1,B1,C1)

def m2DGaussian(x_vals,y_vals,A0,B0,C0,A1,B1,C1):
    return A0*Gauss(x_vals,1,B0,C0)*Gauss(y_vals,1,B0,C0) + A1*Gauss(x_vals,1,B1,C1)*Gauss(y_vals,1,B1,C1)

def FWHM(xs,ys):
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
    a0 = ps[0] / (np.sqrt(2*np.pi) * ps[1]) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = ps[3] / (np.sqrt(2*np.pi) * ps[4]) # Amplitude of second 2D Gaussian (nA/mm^2)
    
    G0,G0_err = integrate.quad(lambda x: Gauss(x,1,ps[1],0),-r,r)
    G1,G1_err = integrate.quad(lambda x: Gauss(x,1,ps[4],0),-r,r)
    
    G2D = a0*G0**2 + a1*G1**2
    
    return G2D

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
    a0 = ps[0] / (np.sqrt(2*np.pi) * ps[1]) # Amplitude of first 2D Gaussian (nA/mm^2)
    a1 = ps[3] / (np.sqrt(2*np.pi) * ps[4]) # Amplitude of second 2D Gaussian (nA/mm^2)
    
    hmpts = np.linspace(-r,r,pts)
    hmxpts,hmypts = np.meshgrid(hmpts,hmpts)
    
    G0x = Gauss(hmxpts,1,ps[1],0)
    G0y = Gauss(hmypts,1,ps[1],0)
    
    G1x = Gauss(hmxpts,1,ps[4],0)
    G1y = Gauss(hmypts,1,ps[4],0)
    
    G = a0*G0x*G0y + a1*G1x*G1y
    
    hm = np.max(G) / 2
    sig2 = ps[4] * 5
    y_sig2 = mGaussianSum(sig2,ps[0],ps[1],ps[2],ps[3],ps[4],ps[5])
    
    contour = plt.contour(hmxpts, hmypts, G, levels=[y_sig2,hm], colors='white', linestyles='dashed', linewidths=2)
    plt.clabel(contour, inline=True, fontsize=8, fmt={y_sig2: '5σ',hm: 'FWHM'})
    
    plt.imshow(G,origin='lower', cmap='viridis',extent=[-r,r,-r,r],vmin=0,vmax=1.5)
    
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')
    plt.title(f'Sum of 2D-Gaussians Fit {mlabel}')
    
    cbar = plt.colorbar()  # Create a colorbar
    cbar.set_label('Current Density (nA/mm^2)', rotation=270, labelpad=15)
    
    if mlabel and save:
        mlabel = '_'.join(mlabel.split())
        plt.savefig(f'./plots/SrSrcs/{mlabel}_HeatMap.png',format='png',dpi=400)
        plt.close()
    else:
        plt.show()
