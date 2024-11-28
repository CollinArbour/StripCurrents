import matplotlib.pyplot as plt
import src.helpers as hp

# Example Usage
log_files = [
    './data/LogFiles/241018_241024_CAENGECO2020.log',
    './data/LogFiles/241105_241112_CAENGECO2020.log',
    './data/LogFiles/241115_241118_CAENGECO2020.log',
    './data/LogFiles/241121_241126_CAENGECO2020.log'
]

start_date = '2024-11-18'
end_date = '2024-11-27'

acc_charge = None
max_charge = 0

for log in log_files:

    matches = hp.parse_log(log)
    timestamps = hp.timestamps(matches)
    imon_values = hp.imon_values(matches)
    
    acc_charge = (hp.accCharge_calc(timestamps, imon_values))
    max_charge += acc_charge[-1]
    #plt.plot(timestamps,acc_charge)

print(max_charge)
plt.show()
