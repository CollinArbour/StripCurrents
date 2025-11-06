import re

def parseRates(txtfl,toGrab):
    with open(txtfl,'r') as fl:
        lines = fl.readlines()
    
    ratePattern = re.compile(r"(\d\w+):")

    rates = {}
    for line in lines:
        if rateLabel := re.search(ratePattern,line):
            if rateLabel[1] in toGrab:
                # Need to incorporate this in regex expression
                rates[rateLabel[1]] = float(line.strip().split()[-1])/10
    
    return rates