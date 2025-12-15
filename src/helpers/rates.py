import re

def parseRates(txtfl,toGrab):
    '''
    Extracts specific rate from a 10s TMB Dump

    Args:
        txtfl: Path to TMB Dump text file
        toGrab: List of rates to extract, who's values shoudl be extracted
    
    Returns:
        Dictionary with rates and their values 
    '''
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