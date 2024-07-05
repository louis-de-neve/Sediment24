import xarray as x
import os
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from operator import itemgetter
import scipy.optimize, scipy.integrate

#IMPORT DATA
luderitz_lat_index = 15
luderitz_long_index = 15

time_array = []
data_array = []

for f in os.listdir("PAR/DATA/requested_files"):
    date = datetime.strptime(f.split(".")[1], '%Y%m%d')
    ds = x.open_dataset(f"PAR/DATA/requested_files/{f}")
    data_point = ds["par"][15, 15].data
    
    time_array.append(date.timetuple().tm_yday)
    data_array.append(data_point)

data_copy = np.asarray(data_array).tolist()
cleanedList = [x for x in data_copy if str(x) != 'nan']
original_mean = np.mean(cleanedList)

time_array, data_array = [list(x) for x in zip(*sorted(zip(time_array, data_array), key = itemgetter(0)))]

#FIRST PASS PAR
def PAR(t):
    return 19*np.sin(0.0172*t + 1.72) + 49


#PLOTTING
fig, axs = plt.subplots(3, 1, figsize=(9, 18), sharex=True)
ax, ax2, ax3 = axs
first_pass_label = r"First Fit: PAR($t$) = $19$ sin$(2 \pi / 365$ $t + 1.72) + 49$"
ax.scatter(time_array, data_array, label="Aqua-MODIS data")
ax.plot(time_array, PAR(np.asarray(time_array)), color="red", label=first_pass_label)
ax.legend(loc="lower left", frameon=False)
#ax.set_xlabel("Days since January 1st 2023")
ax.set_ylabel(r"PAR / einstein m$^{-2}$ day$^{-1}$")
ax.set_ylim(0, 70)
ax.set_xlim(0, 365)
ax.text(290, 2.3, f'Annual Mean = {round(original_mean, 2)}')

# drops erroneous PAR readings
new_time_array = []
new_data_array = []

for i, value in enumerate(time_array):
    if np.abs(data_array[i] - PAR(value)) < 4:
        new_data_array.append(data_array[i])
        new_time_array.append(time_array[i])



        
#optimises
guess = [16, 1.7, 47]

def sinfunc(t, A, p, c):  return A * np.sin(0.0172*t + p) + c
popt, pcov = scipy.optimize.curve_fit(sinfunc, new_time_array, new_data_array, p0=guess)
A, p, c = popt

#REPLOTS WITH NEW DATA
ax2.scatter(new_time_array, new_data_array, label="Aqua-MODIS data")
second_pass_label = r"Second Fit (Optimized Curve): PAR($t$) = $18.9±0.1$ sin$(2 \pi / 365$ $t + 1.711±0.007) + 50.01±0.01$"
ax2.plot(time_array, sinfunc(np.asarray(time_array), A, p, c), color="orange", label=second_pass_label)
ax2.legend(loc="lower left", frameon=False)
#ax2.set_xlabel("Days since January 1st 2023")
ax2.set_ylabel(r"PAR / einstein m$^{-2}$ day$^{-1}$")
ax2.set_ylim(20, 70)
ax2.set_yticks([25, 35, 45, 55, 65])
#ax2.set_xlim(0, 365)

ax3.scatter(time_array, data_array, label="Aqua-MODIS data")
third_pass_label = r"Third Fit (Conserving Annual Mean): PAR($t$) = $18.9±0.1$ sin$(2 \pi / 365$ $t + 1.711±0.007) + 46.62$ "
ax3.plot(time_array, sinfunc(np.asarray(time_array), A, p, original_mean), color="green", label=third_pass_label)
ax3.legend(loc="lower left", frameon=False)
ax3.set_xlabel("Days since January 1st 2023")
ax3.set_ylabel(r"PAR / einstein m$^{-2}$ day$^{-1}$")
ax3.set_ylim(0, 70)
ax3.set_xlim(0, 365)



#PLOTTING STUFF

fig.suptitle("First, Second and Third order fitting of the \nPhotosynthetically Active Radiation in Luderitz, Namibia", fontsize=18)
fig.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig("PAR/ParThirdPass.png")

