# Kiethley Picoammeter Measurements

This script is supposed to be a one stop shop for processing and analyzing data from the Kiethley picoammeter. General uses are to make HV scans (for directly measuring gas gain), or strip scans (for characterizing beam spots, or in general).

In particular the script is meant to be used iteratively.

Modus operandi: A series of measurements are made (your signal), and generally you make some sort of reference measurment to account for you background
1. Signal file
2. Reference file
    + Dark rates, or 0V (?)

## Strip Scan
To do a very simple strip scan at a specific HV point, use this command: 
```
./scripts/mkStripScan.py -f 251205_StripCurrent_HV3600 --hv 3600 --stripScan
```


# To do:
+ [ ] Exapnd to HV scans too!
    + [ ] Change default dataDir (make it dynamic if possible! )