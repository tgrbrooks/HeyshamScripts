import numpy as np
import matplotlib.pyplot as plt
import csv
import math

# Number of free protons in watchman
# Mass [kt] * N free protons in kt of H2O / 10^32
free_prot_per_kT = 0.668559 # [x10^32 free protons/kT H2O]

# Calculate the total watchman mass
watchman_r = 10 # [m]
watchman_h = 2 * watchman_r # [m]
watchman_volume = math.pi * watchman_r**2 * watchman_h # [m^3]
vol_to_mass = 997 * 1e-3 * 1e-3 # 997 kg/m^3 * 1e-3 T/kg * 1e-3 kT/T
watchman_mass = watchman_volume * 1e-3 # [kT]
print("watchman mass = ",watchman_mass," kT")

# Calculate the watchman fiducial mass
watchman_fv_r = 6.7 # [m]
watchman_fv_h = 2 * watchman_fv_r # [m]
watchman_fv_volume = math.pi * watchman_fv_r**2 * watchman_fv_h # [m^3]
watchman_fv_mass = watchman_fv_volume * 1e-3 # [kT]
print("watchman fiducial mass = ",watchman_fv_mass," kT")
fv_ratio = watchman_fv_mass / watchman_mass

# Number of protons in watchman
watchman_p = watchman_mass * free_prot_per_kT # x10^32 free protons
print("watchman free protons = ",watchman_p," x10^32")

# Function to read in rate from geoneutrinos output
def read_rate(file_name, free_protons):
    rate = []
    with open(file_name, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader)
        for row in reader:
            if len(row) < 2:
                continue
            total_tnu = float(row[0]) #interactions/10^32 free protons/year/keV
            watchman_sec_rate = total_tnu * free_protons / 31556952
            rate.append(watchman_sec_rate)
    rate = np.array(rate)
    return rate

# Energy binning
bin_centers = np.linspace(1.805, 9.995, 820) # [MeV]
bin_width = np.mean(np.diff(bin_centers))*1000 # MeV to KeV
print("bin width = ",bin_width, ' KeV')

# Get the different neutrino rates
# Heysham + neutrino background (Hartlepool off)
rate_on = read_rate('reactor_on.csv', watchman_p)
# Hartlepool + neutrino background (Heysham off)
rate_hart_on = read_rate('reactor_hartlepool_onecore_max.csv', watchman_p)
# Neutrino background (Heysham and Hartlepool off)
rate_off = read_rate('reactor_off.csv', watchman_p)

# Get the reactor only rates by subtracting the background
heysham_rate = np.subtract(rate_on, rate_off)
hartlepool_rate = np.subtract(rate_hart_on, rate_off)

# Plot the spectra
plt.plot(bin_centers, heysham_rate, label='heysham')
plt.plot(bin_centers, hartlepool_rate, label='hartlepool')
plt.plot(bin_centers, rate_off, label='$\\nu$ background')
plt.xlabel('Energy [MeV]')
plt.ylabel('Rate [interactions/sec/MeV]')
plt.legend()
plt.show()

# Calculate the integrated rates
total_heysham_rate = sum(heysham_rate)*bin_width
total_hartlepool_rate = sum(hartlepool_rate)*bin_width
total_bkg_rate = sum(rate_off)*bin_width

print('Rates (inner + veto):')
print('Heysham = ',total_heysham_rate,' IBD/s')
print('Hartlepool = ',total_hartlepool_rate,' IBD/s')
print('Background = ',total_bkg_rate,' IBD/s')
print('Rates (inner):')
print('Heysham = ',total_heysham_rate*fv_ratio,' IBD/s')
print('Hartlepool = ',total_hartlepool_rate*fv_ratio,' IBD/s')
print('Background = ',total_bkg_rate*fv_ratio,' IBD/s')

reco_eff = 0.9
R_sig = reco_eff*total_heysham_rate*fv_ratio*60*60*24
other_bkg = (32.4+13.73+13+2.6)/365.25
R_bkg = reco_eff*total_bkg_rate*fv_ratio*60*60*24 + other_bkg

# Discovery sensitivity for uncertain background
# Assuming systematic uncertainty is a flat percentage (x) of background
def significance(dwell, x):
    s = dwell * R_sig
    b = dwell * R_bkg
    return s/np.sqrt(b+np.power(x*b,2))

# Mean dwell time for a given significance (z sigma) and background uncertainty (x)
def dwell_time(z, x):
    return R_bkg/(np.power(R_sig/z, 2)-np.power(x*R_bkg, 2))


# Calculate significances at range of dwell times for different background uncertainties
x = np.linspace(1, 600, 600)
y_a = significance(x, 0)
y_b = significance(x, 0.1)
y_c = significance(x, 0.05)

plt.plot(x, y_a, label='known bkg')
plt.plot(x, y_b, label='10% uncertainty')
plt.plot(x, y_c, label='5% uncertainty')
plt.legend()
plt.xlabel('Dwell time [days]')
plt.ylabel('Significance [N$\sigma$]')
plt.show()

# Calculate dwell times for range of significances for different background uncertainties
z = np.linspace(1, 2.5, 50)
dt_a = dwell_time(z, 0)
dt_b = dwell_time(z, 0.1)
dt_c = dwell_time(z, 0.05)
plt.plot(z, dt_a, label='known bkg')
plt.plot(z, dt_b, label='10% uncertainty')
plt.plot(z, dt_c, label='5% uncertainty')
plt.legend()
plt.xlabel('Significance [N$\sigma$]')
plt.ylabel('Dwell time [days]')
plt.show()
