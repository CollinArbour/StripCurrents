#!/usr/bin/env python3
import re
import numpy as np
import subprocess
from concurrent.futures import ProcessPoolExecutor

def run_macro(args):
    day, hole = args
    hv = 3600
    rebin = 256
    cmd = [
        "root",
        "-l",
        "-b",
        "-q",
        f'fitPeaksGG.C("{day}",{hv},{hole},{rebin})'
        #f'fitSimpleGauss.C("{day}",{hv},{hole},{rebin})'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return (day, hole, result.returncode, result.stdout, result.stderr)

def parsePeak(lines):
    for line in lines:
        if '(p1)' in line:
            return float(line.split()[-1])
    
def getPeakRatios(holes,day,fit,hv,rbf):
    numfl = f'./data_processed/{day}/fitPlots/{fit}/hv{hv}_hole{holes[0]}_rbf{rbf}.txt'
    denfl = f'./data_processed/{day}/fitPlots/{fit}/hv{hv}_hole{holes[1]}_rbf{rbf}.txt'

    num = parsePeak(open(numfl).readlines())
    den = parsePeak(open(denfl).readlines())
    return num/den

if __name__ == "__main__": 
    days = [
        '20240827_dyn_recupCF4',
        # '20240916_dyn_recupCF4', # No source measurements
        '20241025_dyn_recupCF4',
        '20241031_dyn_recupCF4',
        '20241129_dyn_recupCF4',
        # '20241202_dyn_recupCF4', # No source measurements
        '20250121_dyn_recupCF4',
        '20250206_dyn_recupCF4',
        '20250221_dyn_recupCF4',
        '20250312_dyn_recupCF4',
        '20250321_dyn_recupCF4',
        '20250327_dyn_recupCF4',
        '20250401_dyn_recupCF4',
        '20250415_dyn_recupCF4',
        '20250429_dyn_recupCF4',
        '20250505_dyn_recupCF4',
        '20250516_dyn_recupCF4',
        '20250606_dyn_recupCF4',
        '20250617_dyn_recupCF4',
        '20250701_dyn_recupCF4',
        '20250711_dyn_recupCF4',
        '20250729_dyn_recupCF4',
        '20250808_dyn_recupCF4',
        '20250814_dyn_recupCF4',
        '20250821_dyn_recupCF4',
        '20250827_dyn_recupCF4'
    ]
    # Not all holes:
        # '20240916_dyn_recupCF4',
        # '20241202_dyn_recupCF4',
        # '20241121_dyn_recupCF4',
        # '20241122_dyn_recupCF4',

    holes = [2, 4, 6]
    jobs = [(day, hole) for day in days for hole in holes]

    theDate = re.compile(r"(\d{4})(\d{2})(\d{2})")

    with ProcessPoolExecutor() as executor:
        results = executor.map(run_macro, jobs)

    for day, hole, code, out, err in results:
        print(f"🧪 {day}, hole {hole}: Return code {code}")
        if code != 0:
            print(f"❌ stderr:\n{err}")
        # print(out)
    
    fit = 'TwoGauss'
    hv = 3600
    rbf = 256

    # peakPos = []
    # for hole in holes:
    #     peakPos.append([])

    # for day in days:
    #     for i,hole in enumerate(holes):
    #         mfl = f'./data_processed/{day}/fitPlots/{fit}/hv{hv}_hole{hole}_rbf{rbf}.txt'
    #         peakPos[i].append(parsePeak(open(mfl).readlines()))
    #         print(f"🧪 {day}, hole {hole}: {peakPos[i][-1]}")
 
    # peakPos = np.array(peakPos)

    # avgDenom = 0.5 * (peakPos[1] + peakPos[2])
    # ratios = peakPos[0] / avgDenom

    # with open(f'data_processed/ratio_2_46.txt', 'w') as fl:
    #     for val in ratios:
    #         fl.write(f"{val}\n")

    ratios = [(2,4), (2,6)]
    relGasGain = [[] for _ in ratios]


    for day in days:
        for i,ratio in enumerate(ratios):
            relGasGain[i].append(getPeakRatios(ratio, day, fit, hv, rbf))
            print(f"🧪 {day}, holes {ratio}: {relGasGain[i][-1]}")
    
    for i,ratio in enumerate(relGasGain):
        with open(f'data_processed/ratio_{ratios[i][0]}_{ratios[i][1]}.txt', 'w') as fl:
            for j,val in enumerate(ratio):
                fl.write(f"{'-'.join(theDate.search(days[j]).groups())}\t{val}\n")