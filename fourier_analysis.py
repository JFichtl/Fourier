import pandas as pd
import numpy as np
from scipy.fft import rfft, rfftfreq
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

semi_diurnal_species = ['M2', 'S2', 'N2', 'nabla2', 'mu2', '2"N2', 'lambda2', 'T2', 'R2', '2SM2', 'L2', 'K2']
semi_diurnal_periods = np.array([12.4206012, 12, 12.65834751, 12.62600509, 12.8717546, 12.90537297, 12.22177348, 12.01644934, 11.98359564, 11.60695157, 12.19162085,11.96723606])
semi_diurnal_dict = {species : period for species, period in zip(semi_diurnal_species, semi_diurnal_periods)}

diurnal_species = ['K1', 'O1', 'OO1', 'S1', 'M1', 'J1', 'RHO', 'Q1', '2Q1', 'P1']
diurnal_periods = np.array([23.93447213, 25.81933871, 22.30608083, 24, 24.84120241, 23.09848146, 26.72305326, 26.868350, 28.00621204, 24.0658876, 24.06588766])
diurnal_dict = {species : period for species, period in zip(diurnal_species, diurnal_periods)}

long_species = ['Mm', 'Ssa', 'Sa', 'MSr', 'Mr']
long_periods = [661.3111655, 4383.076325, 8766.15265, 354.3670666, 327.8599387]
long_dict = {species : period for species, period in zip(long_species, long_periods)}

short_species = ['M4', 'M6', 'MK3', 'S4', 'MN4', 'S6', 'M3', '2"MK3', 'M8', 'MS4']
short_periods = np.array([6.210300601, 4.140200401, 8.177140247, 6, 6.269173724, 4, 8.280400802, 8.38630265, 3.105150301, 6.103339275])
short_dict = {species : period for species, period in zip(short_species, short_periods)}


def freq_analysis(data: np.ndarray, sample_rate: int, peak_sensitivity: float, unit: str, plot = True):
    N = data.size
    if sample_rate == 0:
        print("sample rate can't be zero.")
        return 0

    xf = rfftfreq(N, 1/sample_rate)
    yf = rfft(data)
    yf[0] = 0

    z = 1.0/N * np.abs(yf)
    h = peak_sensitivity*np.abs(z).max()
    peaks, _ = find_peaks(z, height = h)
    #print(xf[peaks])
    if plot:
        plt.plot(xf, z)
        for k in range(len(peaks)):
            n = peaks[k]
            #plt.scatter(xf[n], z[n], marker = 'x', label = round(xf[n],6))
            plt.scatter(xf[n], z[n], marker = 'x')
        plt.scatter(xf[peaks], z[peaks], marker = 'x')
        plt.xlabel('frecuencia en 1/' + unit)
        plt.ylabel('magnitud')
        #plt.legend()
        plt.show()
    else:
        return xf, z, xf[peaks]

def compare_periods(tide_table, reference):
    measured = return_hourly_freqs(tide_table)
    comparison_dict = {}
    if reference == 'diurnal':
        reference_dict = diurnal_dict
    elif reference == 'semi_diurnal':
        reference_dict = semi_diurnal_dict
    elif reference == 'long':
        reference_dict = long_dict
    elif reference == 'short':
        reference_dict = short_dict
    else:
        print("invalid comparison parameter.")
        return {}

    for species, period in reference_dict.items():
        check = abs(measured - period)
        if len(check) > 0:
            best_index = np.where(check == check.min()) 
            comparison_dict.update({species : measured[best_index][0]})
    return comparison_dict

def import_uhslc_tides(filename):
    tides = pd.read_csv(filename)
    tides.columns = ['year', 'month', 'day', 'hour', 'height']
    mean_height = int(np.nanmean(tides.height))
    tides.height = tides.height.replace(-32767, int(np.nanmean(tides.height)))
    tides['adj_height'] = tides.height - mean_height
    return tides

def return_hourly_freqs(tide_table, threshold = .001, peaks = False):
    freqs, magns, peak_freqs = freq_analysis(tide_table.adj_height, 24, threshold, 'dias', plot = False)
    if peaks:
        return 1/peak_freqs * 24
    else:
        return 1/freqs[1:] * 24

mazatlan_tides = import_uhslc_tides('mazatlan_h.csv')
mazatlan_diurnal_tides = compare_periods(mazatlan_tides, 'diurnal')
