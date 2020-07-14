import numpy as np
import matplotlib.pyplot as plt
import h5py
from bmtk.analyzer.cell_vars import _get_cell_report, plot_report
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.signal import find_peaks
import pickle


#strengths = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 50]
strengths = [0, 2, 4, 6, 8, 10]
f = "percent_results.pkl"
file = open(f, 'rb')
res = pickle.load(file)
file.close()

percents = [res[key] for key in strengths]

plt.figure()
plt.scatter(strengths, percents, marker='.', label='1inh')
plt.title("Exc Synapse Conductance vs Percent AP Conversion")
plt.xticks(strengths)
plt.xlabel("Exc Synapse Condunctance (1e-3)")
plt.ylabel("Percent exc pre-synaptic spikes that caused APs")

strengths = [0, 2, 4, 6, 8, 10]
f = "percent2_results.pkl"
file = open(f, 'rb')
res = pickle.load(file)
file.close()

percents = [res[key] for key in strengths]

#plt.figure()
plt.scatter(strengths, percents, marker='.', label='2inh')
plt.title("Exc Synapse Conductance vs Percent AP Conversion 2s")
plt.xticks(strengths)
plt.xlabel("Exc Synapse Condunctance (1e-3)")
plt.ylabel("Percent exc pre-synaptic spikes that caused APs")
plt.legend()

plt.show()