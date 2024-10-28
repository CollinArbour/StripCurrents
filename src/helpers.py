import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as integrate


import pandas as pd
import re
from datetime import datetime
from scipy.interpolate import interp1d


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


def parse_log(log_data):
    pattern = r"\[(.*?)\]:.*par \[IMonH\] val \[(.*?)\];"
    matches = re.findall(pattern, log_data)
    
    return matches

def timestamps(matches):
    timestamps = [datetime.fromisoformat(match[0]) for match in matches]

    return timestamps

def imon_values(matches):
    imon_values = [float(match[1]) for match in matches]
    
    return imon_values

def current_vs_time(start_date, end_date, timestamps, imon_values):
    
    df = pd.DataFrame({'Timestamp': timestamps, 'IMon': imon_values})
    df.set_index('Timestamp', inplace=True)

    # Step 3: Filter the data based on the provided date range
    df_filtered = df[start_date:end_date]

    # Step 4: Resample data to reduce noise (e.g., take mean every minute)
    df_resampled = df_filtered.resample('5min').mean()  # Resample per minute, adjust '1T' for different intervals (e.g., '5T' for 5 minutes)

    #print(df_resampled.mean()* 0.7/1.2)

    #dfData = (df_resampled)
    #imon = open('./imonvalues.txt', 'a')
    #imon.write(df_resampled.to_string())

    #imon = pd.DataFrame.to_string
    #imon.write(imon)


    #print(df_resampled)

    # Step 5: Plot the resampled data
    plt.plot(df_resampled.index, df_resampled['IMon'], label="IMon (μA)", marker=".", color="blue")

    max_imonh = df_resampled['IMon'].max()  # Find the maximum IMonH value
    top_limit = max_imonh * 1.2  # Increase the y-axis limit by 20%


    # Customize the plot
    plt.title(f"IMon Over Time ({start_date} to {end_date})")
    plt.xlabel("Time (Days)")
    plt.ylabel("IMon (μA)")
    #plt.grid(True)
    plt.xticks(rotation=45)
    plt.ylim(0, top_limit)
    plt.tight_layout()
    plt.legend()

    # Show plot
    plt.show()

# This is the skeleton of how we will be getting our graph for CAEN current vs time
'''
def current_vs_time(start_date, end_date):
    # Step 1: Read the data from a text file
    file_path = './data/LogFiles/CAENGECO2020.log'

    with open(file_path, 'r') as file:
        data = file.read()

    # Step 2: Parse data for timestamps and IMonH values
    pattern = r"\[(.*?)\]:.*par \[IMonH\] val \[(.*?)\];"
    matches = re.findall(pattern, data)

    # Convert the extracted data into a DataFrame
    timestamps = [datetime.fromisoformat(match[0]) for match in matches]
    imon_values = [float(match[1]) for match in matches]

    df = pd.DataFrame({'Timestamp': timestamps, 'IMon': imon_values})
    df.set_index('Timestamp', inplace=True)

    # Step 3: Filter the data based on the provided date range
    df_filtered = df[start_date:end_date]

    # Step 4: Resample data to reduce noise (e.g., take mean every minute)
    df_resampled = df_filtered.resample('5min').mean()  # Resample per minute, adjust '1T' for different intervals (e.g., '5T' for 5 minutes)

    print(df_resampled.mean()* 0.7/1.2)

    #dfData = (df_resampled)
    #imon = open('./imonvalues.txt', 'a')
    #imon.write(df_resampled.to_string())

    #imon = pd.DataFrame.to_string
    #imon.write(imon)


    #print(df_resampled)

    # Step 5: Plot the resampled data
    plt.plot(df_resampled.index, df_resampled['IMon'], label="IMon (μA)", marker=".", color="blue")

    max_imonh = df_resampled['IMon'].max()  # Find the maximum IMonH value
    top_limit = max_imonh * 1.2  # Increase the y-axis limit by 20%


    # Customize the plot
    plt.title(f"IMon Over Time ({start_date} to {end_date})")
    plt.xlabel("Time (Days)")
    plt.ylabel("IMon (μA)")
    #plt.grid(True)
    plt.xticks(rotation=45)
    plt.ylim(0, top_limit)
    plt.tight_layout()
    plt.legend()

    # Show plot
    plt.show()
'''

def accCharge_vs_time():
    print()

    log_file_path ='./data/LogFiles/CAENGECO2020.log'

    with open(log_file_path, 'r') as file:
        log_data = file.read()

    pattern = r"\[(.*?)\]:.*par \[IMonH\] val \[(.*?)\];"
    matches = re.findall(pattern, log_data)

    timestamps = [datetime.fromisoformat(match[0]) for match in matches]
    imon_values = [float(match[1]) for match in matches]



# instead of I vs T we are going to do accQ vs T
def accCharge_vs_time(start_date, end_date):

    #take in accumulated charge as a value and plot it versus time 

    start_date = datetime.strptime(start_date, '%Y-%m-%d').date()
    end_date = datetime.strptime(end_date, '%Y-%m-%d').date()

    log_file_path = './data/LogFiles/CAENGECO2020.log'
    acc_charge_file_path = './data/accChrg_vtime.txt'

    # Step 1: Read log file and parse full timestamps
    with open(log_file_path, 'r') as file:
        log_data = file.read()

    # Regex to extract full timestamps and IMonH values
    pattern = r"\[(.*?)\]:.*par \[IMonH\] val \[(.*?)\];"
    matches = re.findall(pattern, log_data)

    # Convert matches to list of datetime objects for timestamps
    timestamps = [datetime.fromisoformat(match[0]) for match in matches]

    # Step 2: Load accumulated charge values
    with open(acc_charge_file_path, 'r') as file:
        acc_charge_values = [float(line.strip()) for line in file if line.strip()]

    #print(acc_charge_values)

    # Interpolate the accumulated charge values to match the number of timestamps
    original_indices = np.linspace(0, len(acc_charge_values) - 1, num=len(acc_charge_values))
    target_indices = np.linspace(0, len(acc_charge_values) - 1, num=len(timestamps))

    interpolation_function = interp1d(original_indices, acc_charge_values, kind='linear')
    resampled_acc_charge_values = interpolation_function(target_indices)

    # Create DataFrame with matched timestamps and resampled accumulated charge values
    df = pd.DataFrame({'Timestamp': timestamps, 'Accumulated Charge': resampled_acc_charge_values})
    df.set_index('Timestamp', inplace=True)

    # Step 3: Filter data based on start_date and end_date
    df_filtered = df.loc[(df.index.date >= start_date) & (df.index.date <= end_date)]

    max_acc_charge = df_filtered['Accumulated Charge'].max()
    top_limit = max_acc_charge * 1.2

    # Plot precise data
    plt.plot(df_filtered.index, df_filtered['Accumulated Charge'], label="Accumulated Charge", color="blue")
    plt.title(f"Accumulated Charge Over Time ({start_date} to {end_date})")
    plt.text(0.25, 0.94, f'Max Accumulated Charge Reached: {max_acc_charge:.2f}',
         horizontalalignment='center',
         verticalalignment='center',
         transform=plt.gca().transAxes,
         fontsize=8,
         bbox=dict(facecolor='white', alpha=0.5))
    plt.xlabel("Time (Days)")
    plt.ylabel("Accumulated Charge (mC/cm)")
    plt.xticks(rotation=45)
    plt.ylim(0, top_limit)
    plt.xlim(left=start_date)
    plt.tight_layout()
    plt.legend()
    plt.show()


