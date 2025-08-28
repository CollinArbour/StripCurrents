#!/usr/bin/env python3
import re
import numpy as np
import subprocess
from concurrent.futures import ProcessPoolExecutor

from load_config import config

def run_macro(args):
    day, hole = args
    hv = 3600
    rebin = 256
    cmd = [
        "root",
        "-l",
        "-b",
        "-q",
        f'src/fitPeaksGG.C("{day}",{hv},{hole},{rebin})'
        #f'fitSimpleGauss.C("{day}",{hv},{hole},{rebin})'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return (day, hole, result.returncode, result.stdout, result.stderr)

def parsePeak(lines):
    for line in lines:
        if '(p1)' in line:
            return float(line.split()[-1])
    
def getPeakRatios(holes,day,fit,hv,rbf):
    numfl = f'data/GasGain/data_processed/{day}/fitPlots/{fit}/hv{hv}_hole{holes[0]}_rbf{rbf}.txt'
    denfl = f'data/GasGain/data_processed/{day}/fitPlots/{fit}/hv{hv}_hole{holes[1]}_rbf{rbf}.txt'

    num = parsePeak(open(numfl).readlines())
    den = parsePeak(open(denfl).readlines())
    return num/den

if __name__ == "__main__": 
    holes = [2, 4, 6]
    jobs = [(day, hole) for day in config.measurementDays for hole in holes]

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

    ratios = [(2,4), (2,6)]
    relGasGain = [[] for _ in ratios]

    for day in config.measurementDays:
        for i,ratio in enumerate(ratios):
            relGasGain[i].append(getPeakRatios(ratio, day, fit, hv, rbf))
            print(f"🧪 {day}, holes {ratio}: {relGasGain[i][-1]}")
    
    for i,ratio in enumerate(relGasGain):
        with open(f'data/GasGain/data_processed/ratio_{ratios[i][0]}_{ratios[i][1]}.txt', 'w') as fl:
            for j,val in enumerate(ratio):
                fl.write(f"{'-'.join(theDate.search(config.measurementDays[j]).groups())}\t{val}\n")