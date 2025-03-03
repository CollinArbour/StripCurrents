#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

dates = ['20240827',
        '20241025',
        '20241031',
        '20241121',
        '20241129',
        '20250121',
        '20250206']

# Hole2
#days = ['20240827_dyn_recup_CF4',
#        '20241025_dyn_recupCF4',
#        '20241031_dyn_recupCF4',
#        '20241121_dyn_recupCF4',
#        '20241129_dyn_recupCF4',
#        '06_dyn_recupCF4',
#        '07_dyn_recupCF4']


# Hole4 and Hole6
days = ['20240827_dyn_recup_CF4',
        '20241025_dyn_recupCF4',
        '20241031_dyn_recupCF4',
        '20241122_dyn_recupCF4',
        '20241129_dyn_recupCF4',
        '06_dyn_recupCF4',
        '07_dyn_recupCF4']

#for hole in [1,3,5]:
#    print(f'for hole: {hole}')
#    for day in days:
#        path = f'data_processed.nosync/{day}/wSrc/hole{hole}/fitParams.txt'
#
#        try:
#            with open(path,'r') as fl:
#                lines = fl.readlines()
#        except:
#            print(f'\tFile not found: {path}')
#            continue
#        
#        date = []
#        pos = []
#        
#        date.append(day)
#        for line in lines:
#            if 'Position' in line:
#                pos.append(float(line.split(':')[-1]))
#
#        print(f'\t{day},{pos}')


accChrg= np.array([0.0,11.68,11.68,31.13,48.38,61.68,98.02])

##Optained from previous run iteration
### Hole 2 --------------------------------------------------
peakPos_hole2 = [1744.0897098607102,1762.934905544564,1543.5746213907626,2013.0461815805709,1684.9920610511986,1941.6802834613113,1693.7676685442764]
#plt.annotate('966.8\n mbar',(accChrg[0],peakPos_hole2[0]),textcoords='offset points',xytext=(0,-30),ha='center')
#plt.annotate('962.1 mbar',(accChrg[1],peakPos_hole2[1]),textcoords='offset points',xytext=(0,10),ha='center')
#plt.annotate('970.6 mbar',(accChrg[2],peakPos_hole2[2]),textcoords='offset points',xytext=(40,-5),ha='center')
#plt.annotate('944.8 mbar',(accChrg[3],peakPos_hole2[3]),textcoords='offset points',xytext=(40,-5),ha='center')
#plt.annotate('955.5 mbar',(accChrg[4],peakPos_hole2[4]),textcoords='offset points',xytext=(0,-30),ha='center')
#
### Other holes --------------------------------------------------
offMask = np.array([True,False,True,False,False,True,True])
peakPos_hole1 = np.array([1719.700766522691,1514.3171954717113,1887.2633419566137,1683.1422288483657])
peakPos_hole3 = np.array([1673.9549701047274,1506.129806205713,1912.9893997745803,1682.1798984063726])
peakPos_hole5 = np.array([1713.4261548300008,1560.8827974687767,1916.6469563875391,1733.6706102880507])

peakPos_hole4 = [1741.2082639872517,1758.9679977006772,1523.8945529231733,1949.3193016609603,1657.0404169135736,1901.6861399748645,1702.838069886899]
peakPos_hole6 = [1603.9568116270716,1684.4408097883668,1434.5605298881815,1836.8142504691432,1579.9764335042576,1837.0636044001221,1642.212066780308]

### Plot Ratios --------------------------------------------------
ratio24 = np.array(peakPos_hole2)/np.array(peakPos_hole4)
ratio26 = np.array(peakPos_hole2)/np.array(peakPos_hole6)
ratio21 = np.array(peakPos_hole2)[offMask]/np.array(peakPos_hole1)
ratio23 = np.array(peakPos_hole2)[offMask]/np.array(peakPos_hole3)
ratio25 = np.array(peakPos_hole2)[offMask]/np.array(peakPos_hole5)

### Same strip ratio plots
ratio36 = np.array(peakPos_hole3)/np.array(peakPos_hole6)[offMask]
ratio14 = np.array(peakPos_hole1)/np.array(peakPos_hole4)[offMask]
ratio25 = np.array(peakPos_hole2)[offMask]/np.array(peakPos_hole5)

### Plot ratio of peak positions
#plt.plot(accChrg[offMask],ratio21,marker='o',linestyle='',color='tab:blue',label='Hole2/Hole1')
#plt.plot(accChrg[offMask],ratio23,marker='o',linestyle='',color='tab:brown',label='Hole2/Hole3')
#plt.plot(accChrg,ratio24,marker='o',linestyle='',color='tab:green',label='Hole2/Hole4')
##plt.plot(accChrg[offMask],ratio25,marker='o',linestyle='',color='tab:purple',label='Hole2/Hole5')
#plt.plot(accChrg,ratio26,marker='o',linestyle='',color='tab:orange',label='Hole2/Hole6')

plt.plot(accChrg[offMask],ratio36,marker='o',linestyle='',color='tab:blue',label='Hole3/Hole6')
plt.plot(accChrg[offMask],ratio14,marker='o',linestyle='',color='tab:green',label='Hole1/Hole4')
plt.plot(accChrg[offMask],ratio25,marker='o',linestyle='',color='black',label='Hole2/Hole5')


# Raw plots
#plt.plot(accChrg,peakPos_hole2,marker='o',linestyle='',color='black',label='Hole2')
#
#plt.plot(accChrg[offMask],peakPos_hole1,marker='o',linestyle='',color='tab:blue',label='Hole1')
#plt.plot(accChrg[offMask],peakPos_hole3,marker='o',linestyle='',color='tab:brown',label='Hole3')
#plt.plot(accChrg,peakPos_hole4,marker='o',linestyle='',color='tab:green',label='Hole4')
#plt.plot(accChrg[offMask],peakPos_hole5,marker='o',linestyle='',color='tab:purple',label='Hole5')
#plt.plot(accChrg,peakPos_hole6,marker='o',linestyle='',color='tab:cyan',label='Hole6')


plt.legend()
plt.axhline(1,linestyle='--',color='black',alpha=0.3)
plt.xlabel('Accumulated Charge (mC/cm)')

plt.ylabel('Ratio of Peak Positions (ADC)')
plt.title('Ratio of 109Cd Peak Positions')

#plt.ylabel('Peak Position (ADC)')
#plt.title('109Cd Peak Positions')


#plt.savefig('./hole4_peakPos_vs_accChrg.png',format='png',dpi=200)
plt.savefig('./sameStripHoles_peakPosRatio_vs_accChrg.png',format='png',dpi=200)
plt.close()