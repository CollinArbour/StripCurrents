# Strip Current Measurments

## Format for data header:
```
MiniCSC_4
L_1:S_3:HV_0:Imon_0.0
src_No:NoSrc
FlowRate_0.5
GasComp:0.4_Ar:0.55_CO2:0.05_CF4_recup
Comments:
```

## Plot title format:
```
MiniCSC4 {Measurement type} L1 H2 90Sr-Src{i+1} HV3600
        strip scan vs hvSCAN 
```

## To do items:
+ [ ] Split up `helpers.py`
+ [ ] Create data run object for GasGain DAQ  Measurements
  + [ ] Functions for and storing fit data
  + [ ] Save (read and write to text) so it can pick back up