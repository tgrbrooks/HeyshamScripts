import ROOT
import numpy as np
import matplotlib.pyplot as plt
import random
from math import sqrt
from scipy import optimize

# Script for doing a simplified Heysham analysis
# Useful for testing analysis techniques (cut optimisation, energy spectrum 
# fits, etc) outside of watchmakers

# Global variables
num_samples = 100000 # Number of random numbers to generate at various points
detector_r = 10.02 # Detector radius [m]
pmt_r = 6.94 # PMT support structure radius [m]
signal_eff = 0.9 # Signal reconstruction efficiency (I assume due to proximity cut??)

# Interaction rates [per second]
rates = {"IBDPositron_LIQUID_ibd_p" : 9.663984059177409e-05, 
         "IBDNeutron_LIQUID_ibd_n" : 9.663984059177409e-05,
         "IBD_LIQUID_pn_ibd" : 9.663984059177409e-05,
         "heysham_signal" : 9.02e-6,
         "IBDPositronHeyshamSig_LIQUID_ibd_p_hs" : 9.02e-6,
         "heysham_background" : 2.484e-5,
         "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb" : 2.484e-5,
         "228Ac_LIQUID_CHAIN_232Th_NA" : 8.82E-01, #
         "212Pb_LIQUID_CHAIN_232Th_NA" : 8.82E-01, #
         "212Bi_LIQUID_CHAIN_232Th_NA" : 5.65E-01, #
         "208Tl_LIQUID_CHAIN_232Th_NA" : 3.17E-01, #
         "214Pb_LIQUID_CHAIN_222Rn_NA" : 6.3, #
         "214Bi_LIQUID_CHAIN_222Rn_NA" : 6.3,#
         "210Bi_LIQUID_CHAIN_222Rn_NA" : 6.3,#
         "210Tl_LIQUID_CHAIN_222Rn_NA" : 1.26e-3,#
         "40K_LIQUID_40K_NA" : 25.2,#
         "n17_LIQUID_A_Z" : 8.01e-06,
         "li9_LIQUID_A_Z" : 7.97e-06,
         "234Pa_PMT_CHAIN_238U_NA" : 2.49E+03,#
         "214Pb_PMT_CHAIN_238U_NA" : 2.49E+03,#
         "214Bi_PMT_CHAIN_238U_NA" : 2.49E+03,#
         "210Bi_PMT_CHAIN_238U_NA" : 2.49E+03,#
         "210Tl_PMT_CHAIN_238U_NA" : 4.98E-01,#
         "228Ac_PMT_CHAIN_232Th_NA" : 2.54E+03,#
         "212Pb_PMT_CHAIN_232Th_NA" : 2.54E+03,#
         "212Bi_PMT_CHAIN_232Th_NA" : 1.62E+03,#
         "208Tl_PMT_CHAIN_232Th_NA" : 9.13E+02, #
         "40K_PMT_40K_NA" : 5.08E+03} #

# Save keys for everything that isn't an accidental
p_key = "IBDPositron_LIQUID_ibd_p"
n_key = "IBDNeutron_LIQUID_ibd_n"
ibd_key = "IBD_LIQUID_pn_ibd"
hs_key = "IBDPositronHeyshamSig_LIQUID_ibd_p_hs"
hb_key = "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb"
not_acc = [p_key, n_key, ibd_key, hs_key, hb_key]

# Class for storing different interactions 
class SingleEvent:
    def __init__(self, name, rate, evts):
        self.name = name
        self.rate = rate
        self.labels = {"mc_x" : "$x^{MC}$ [m]",
                       "mc_y" : "$y^{MC}$ [m]",
                       "mc_z" : "$y^{MC}$ [m]",
                       "reco_x" : "$x^{reco}$ [m]",
                       "reco_y" : "$y^{reco}$ [m]",
                       "reco_z" : "$z^{reco}$ [m]",
                       "mc_energy" : "$E^{MC}$ [GeV]",
                       "inner_hits" : "Inner PMT hits",
                       "veto_hits" : "Veto PMT hits",
                       "n9" : "n9",
                       "dwall" : "Distance to wall (MC) [m]",
                       "reco_dwall" : "Distance to wall [m]"}
        self.values = {}
        self.values['mc_x'] = []
        self.values['mc_y'] = []
        self.values['mc_z'] = []
        self.values['reco_x'] = []
        self.values['reco_y'] = []
        self.values['reco_z'] = []
        self.values['mc_energy'] = []
        self.values['inner_hits'] = []
        self.values['veto_hits'] = []
        self.values['dwall'] = []
        self.values['reco_dwall'] = []
        self.values['n9'] = []
        self.values['good_pos'] = []
        self.values['closest_pmt'] = []
        self.nevents = evts

    # Add an event to the lists
    def add(self, mc_x, mc_y, mc_z, reco_x, reco_y, reco_z, mc_energy,
            inner_hits, veto_hits, n9, good_pos, closest_pmt):
        self.values['mc_x'].append(mc_x)
        self.values['mc_y'].append(mc_y)
        self.values['mc_z'].append(mc_z)
        self.values['reco_x'].append(reco_x)
        self.values['reco_y'].append(reco_y)
        self.values['reco_z'].append(reco_z)
        self.values['mc_energy'].append(mc_energy)
        self.values['inner_hits'].append(inner_hits)
        self.values['veto_hits'].append(veto_hits)
        self.values['n9'].append(n9)
        self.values['good_pos'].append(good_pos)
        self.values['closest_pmt'].append(closest_pmt)
        # True distance to the PMT support structure
        mc_r = sqrt(mc_x**2 + mc_y**2)
        dwall = 6.7 - max(mc_r, abs(mc_z))
        self.values['dwall'].append(dwall)
        # Reconstructed distance to the PMT support structure
        reco_r = sqrt(reco_x**2 + reco_y**2)
        reco_dwall = 6.7 - max(reco_r, abs(reco_z))
        self.values['reco_dwall'].append(reco_dwall)

    # Replace all with numpy arrays
    def numpify(self):
        for key in self.values:
            self.values[key] = np.array(self.values[key])

    # Reconstruction efficiency
    def efficiency(self, wall_cut, energy_cut, reco=True):
        # with perfect position reconstruction
        eff = np.count_nonzero((self.values['dwall'] > wall_cut)
                               & (self.values['inner_hits'] > energy_cut)
                               & (self.values['veto_hits'] < 4))
        # with actual position reconstruction
        if reco:
            eff = np.count_nonzero((self.values['closest_pmt'] > wall_cut)
                                    & (self.values['n9'] > energy_cut)
                                    & (self.values['good_pos'] > 0.1)
                                    & (self.values['inner_hits'] > 4)
                                    & (self.values['veto_hits'] < 4))
        eff = eff/self.nevents
        return eff

    def eff_rate(self, wall_cut, energy_cut, reco=True):
        return self.rate * self.efficiency(wall_cut, energy_cut, reco)

    # Check if a specific event is selected
    def selected(self, wall_cut, energy_cut, index, reco=True):
        sel = False
        if index >= len(self.values['dwall']):
            return sel
        if not reco and ((self.values['dwall'][index] > wall_cut)
                         & (self.values['inner_hits'][index] > energy_cut)
                         & (self.values['veto_hits'][index] < 4)):
            sel = True
        elif reco and ((self.values['closest_pmt'][index] > wall_cut)
                       & (self.values['n9'][index] > energy_cut)
                       & (self.values['good_pos'][index] > 0.1)
                       & (self.values['inner_hits'][index] > 4)
                       & (self.values['veto_hits'][index] < 4)):
            sel = True
        return sel


# Calculate the reconstructed IBD rate
def calc_signal_rate(singles, wall_p_cut, energy_p_cut, wall_n_cut,
                     energy_n_cut, dt_cut, ds_cut, pos_key, reco=True, 
                     verbose=False, gd=False):

    # Ratio of fiducial volume to detector volume for scaling rates
    fiducial_r = pmt_r - wall_p_cut
    volume_ratio = pow(fiducial_r,3)/pow(detector_r,3)

    # Calculate positron and neutron singles rates and efficiencies
    p_rate = singles[pos_key].eff_rate(wall_p_cut, energy_p_cut, reco)
    n_rate = singles[n_key].eff_rate(wall_n_cut, energy_n_cut, reco)
    p_eff = singles[pos_key].efficiency(wall_p_cut, energy_p_cut, reco)/volume_ratio
    n_eff = singles[n_key].efficiency(wall_n_cut, energy_n_cut, reco)/volume_ratio

    # Draw random time and distance separations from pair distribution
    # Assumption that distance-time distribution not correlated with p,n
    # vertices or energies
    # TODO using fairly simple models for neutron interaction time and distance
    # would be better to just draw from a high stats histogram
    exp_constant = 200
    chisq_val = 8.5
    if gd:
        exp_constant = 25
        chisq_val = 5.5
    sep_time = np.random.exponential(exp_constant, size=num_samples) # [mus]
    # Dist correlated to time, but with perfect position resolution distance
    # cut is not going to remove any pairs so doesn't matter
    sep_dist = np.random.chisquare(chisq_val, size=num_samples) # [cm]

    # Calculate pair rate from time and distance cut
    # TODO not accounting for correlations between time and distance (but the
    # distance is effectively washed out by the reconstruction performance
    # anyway - and so are any correlations)
    time_eff = np.count_nonzero(sep_time < dt_cut)/num_samples
    dist_eff = np.count_nonzero(sep_dist < ds_cut)/num_samples
    time_eff = 1.
    dist_eff = 1.

    # Calculate the total IBD rate [per day]
    signal_rate = p_rate * n_eff * time_eff * dist_eff * 86400
    if verbose:
        print("Positron efficiency = ",p_eff)
        print("Positron rate = ",p_rate*86400,' per day')
        print('Neutron efficiency = ',n_eff)
        print("Neutron rate = ",n_rate*86400,' per day')
        print('Time cut efficiency = ',time_eff)
        print('Distance cut efficiency = ',dist_eff)
        print('Signal rate = ',signal_rate,' per day\n')
    return signal_rate

# Calculate the rate of accidentals
def calc_accidental_rate(singles, wall_p_cut, energy_p_cut, wall_n_cut,
                         energy_n_cut, dt_cut, ds_cut, reco=True,
                         verbose=False):
    # Individual accidental rates
    accidental_p_rates = {}
    accidental_n_rates = {}
    # Loop over all decays
    for key in singles:
        if not ('238U' in key or '232Th' in key or '40K' in key or '222Rn' in key): 
            continue
        # Calculate singles rates from efficiencies for positron and neutron cuts
        accidental_p_rates[key] = singles[key].eff_rate(wall_p_cut,
                                                        energy_p_cut, reco)
        accidental_n_rates[key] = singles[key].eff_rate(wall_n_cut,
                                                        energy_n_cut, reco)
        if verbose:
            print(key,' positron rate = ',accidental_p_rates[key],' per day',
                  ' neutron rate = ',accidental_n_rates[key], ' per day')

    # Sum for total rates
    accidental_p_rate = sum(accidental_p_rates.values())
    accidental_n_rate = sum(accidental_n_rates.values())
    if verbose:
        print('Accidental positron rate = ',accidental_p_rate,' per day',
              ' neutron rate = ',accidental_n_rate,' per day')

    # Calculate distance cut efficiency by drawing a random pair of decays
    # according to reconstructed rates rates
    p_probs = [x/accidental_p_rate for x in list(accidental_p_rates.values())]
    n_probs = [x/accidental_n_rate for x in list(accidental_n_rates.values())]
    rand_p_keys = np.random.choice(list(accidental_p_rates), num_samples,
                                   p=p_probs)
    rand_n_keys = np.random.choice(list(accidental_n_rates), num_samples,
                                   p=n_probs)
    accidental_dist_eff = 0.
    # Loop over the randomly sampled decays
    for i, rp_key in enumerate(rand_p_keys):
        rn_key = rand_n_keys[i]
        # Choose a index at random for each decay
        p_idx = np.random.randint(0, len(singles[rp_key].values['mc_x']), 1)[0]
        # Keep going until a valid interaction is found
        while singles[rp_key].selected(wall_p_cut, energy_p_cut, p_idx, reco):
            p_idx = np.random.randint(0, len(singles[rp_key].values['mc_x']), 1)[0]

        # Do the same for neutron cuts
        n_idx = np.random.randint(0, len(singles[rn_key].values['mc_x']), 1)[0]
        while singles[rn_key].selected(wall_n_cut, energy_n_cut, n_idx, reco):
            n_idx = np.random.randint(0, len(singles[rn_key].values['mc_x']), 1)[0]

        # Get the vertices
        vertex_str = ['mc_x', 'mc_y', 'mc_z']
        if reco:
            vertex_str = ['reco_x', 'reco_y', 'reco_z']
        p_vtx = np.array((singles[rp_key].values[vertex_str[0]][p_idx],
                          singles[rp_key].values[vertex_str[1]][p_idx],
                          singles[rp_key].values[vertex_str[2]][p_idx]))
        n_vtx = np.array((singles[rn_key].values[vertex_str[0]][n_idx],
                          singles[rn_key].values[vertex_str[1]][n_idx],
                          singles[rn_key].values[vertex_str[2]][n_idx]))
        # If the distance between vertices is within the distance cut select
        # the event
        dist = np.linalg.norm(p_vtx-n_vtx)
        if dist < ds_cut/100:
            accidental_dist_eff += 1.
    accidental_dist_eff /= num_samples
    accidental_dist_eff = 0.05

    # Calculate total accidental rate by scaling by time cut and distance
    # efficiency
    accidental_rate = (accidental_p_rate * accidental_n_rate * dt_cut/1e6
                       * accidental_dist_eff * 86400)
    if verbose:
        print('Accidental time cut efficiency = ',dt_cut/1e6)
        print('Accidental distance cut efficiency = ',accidental_dist_eff)
        print('Accidental background rate = ',accidental_rate,' per day\n')
    return accidental_rate


# Calculate the sensitivity in number of sigma (assuming no systematic
# background uncertianty)
def sensitivity(singles, wall_p_cut, energy_p_cut, wall_n_cut, energy_n_cut,
                dt_cut, ds_cut, reco=True, verbose=False, gd=False):
    # Calculate the selected Heysham IBD rate
    signal_rate = calc_signal_rate(singles, wall_p_cut, energy_p_cut,
                                   wall_n_cut, energy_n_cut, dt_cut, ds_cut,
                                   hs_key, reco, verbose, gd)
    # Calculate the accidental rate
    accidental_rate = calc_accidental_rate(singles, wall_p_cut, energy_p_cut,
                                           wall_n_cut, energy_n_cut, dt_cut,
                                           ds_cut, reco, verbose)
    # Calculate the IBD background rate
    ibd_bkg_rate = calc_signal_rate(singles, wall_p_cut, energy_p_cut,
                                    wall_n_cut, energy_n_cut, dt_cut, ds_cut,
                                    hs_key, reco, verbose, gd)
    return signal_rate/sqrt(signal_rate+accidental_rate)


# Calculate the dwell time required to reach a median z sigma detection of
# Heysham by scaling Hartlepool signal (assuming no systematic background
# uncertainty)
def dwell_time_scaling(singles, wall_p_cut, energy_p_cut, wall_n_cut,
                      energy_n_cut, dt_cut, ds_cut, z, reco=True,
                      verbose=False, gd=False):
    # Calculate the selected Hartlepool IBD rate
    signal_rate = calc_signal_rate(singles, wall_p_cut, energy_p_cut,
                                   wall_n_cut, energy_n_cut, dt_cut, ds_cut,
                                   p_key, reco, verbose, gd)
    # Calculate the accidental rate
    accidental_rate = calc_accidental_rate(singles, wall_p_cut, energy_p_cut,
                                           wall_n_cut, energy_n_cut, dt_cut,
                                           ds_cut, reco, verbose)

    # Scale Hartlepool to expected IBD background and add to accidental rate
    background_rate = signal_rate * rates[hb_key]/rates[p_key] * signal_eff + accidental_rate
    # Scale Hartlepool to expected Heysham signal
    signal_rate = signal_rate * rates[hs_key]/rates[p_key] * signal_eff
    if verbose:
        print('Total signal rate = ',signal_rate)
        print('Total background rate = ',background_rate)

    dwell_time = background_rate/pow(signal_rate/z, 2)
    return dwell_time

# Calculate the dwell time required to reach a median z sigma detection of
# Heysham using correct energy spectra (assuming no systematic background
# uncertainty)
def dwell_time_heysham(singles, wall_p_cut, energy_p_cut, wall_n_cut,
                       energy_n_cut, dt_cut, ds_cut, z, reco=True,
                       verbose=False, gd=False):
    # Calculate the selected Heysham IBD rate
    signal_rate = calc_signal_rate(singles, wall_p_cut, energy_p_cut,
                                   wall_n_cut, energy_n_cut, dt_cut, ds_cut,
                                   hs_key, reco, verbose, gd)
    # Calculate the IBD background rate
    ibd_bkg_rate = calc_signal_rate(singles, wall_p_cut, energy_p_cut,
                                    wall_n_cut, energy_n_cut, dt_cut, ds_cut,
                                    hb_key, reco, verbose, gd)
    # Calculate the accidental rate
    accidental_rate = calc_accidental_rate(singles, wall_p_cut, energy_p_cut,
                                           wall_n_cut, energy_n_cut, dt_cut,
                                           ds_cut, reco, verbose)

    signal_rate = signal_rate * signal_eff
    background_rate = ibd_bkg_rate * signal_eff + accidental_rate
    if verbose:
        print('Total signal rate = ',signal_rate)
        print('Total background rate = ',background_rate)

    dwell_time = background_rate/pow(signal_rate/z, 2)
    return dwell_time

# Plot 1D distributions on same canvas
def plot_all(singles, val_key, minmax):
    plt.cla()
    for key in singles:
        plt.hist(singles[key].values[val_key],
                 histtype='step',
                 range=minmax,
                 density=True,
                 bins=100,
                 label=singles[key].name)
        plt.xlabel(singles[key].labels[val_key])
    plt.ylabel('Fraction of events')
    plt.legend(fontsize='xx-small', ncol=2)
    plt.savefig('plots/'+val_key)

root_file = ROOT.TFile.Open("gd_water_reco.root")

# Read in interactions from file
singles = {}
for event in root_file.tree:
    key = str(event.int_name)
    if key not in singles:
        print(key)
        singles[key] = SingleEvent(key, rates[key], event.nevents)
    singles[key].add(event.mc_x/1000., event.mc_y/1000., event.mc_z/1000.,
                     event.reco_x/1000., event.reco_y/1000., event.reco_z/1000.,
                     event.mc_energy, event.inner_hits, event.veto_hits,
                     event.n9, event.good_pos, event.closest_pmt/1000.)


for key in singles:
    singles[key].numpify()

plot_all(singles, 'mc_x', (-10, 10))
plot_all(singles, 'mc_y', (-10, 10))
plot_all(singles, 'mc_z', (-10, 10))
plot_all(singles, 'mc_energy', (0, 10))
plot_all(singles, 'inner_hits', (0, 100))
plot_all(singles, 'veto_hits', (0, 100))
plot_all(singles, 'dwall', (-5, 10))

# Playing around with optimisation
'''
def dwell_time_opt(x):
    signal_rate = calc_signal_rate(singles, x[0], x[1], x[2], x[3], x[4], x[5],
                                   hs_key)
    ibd_bkg_rate = calc_signal_rate(singles, x[0], x[1], x[2], x[3], x[4], x[5],
                                    hb_key)
    accidental_rate = calc_accidental_rate(singles, x[0], x[1], x[2], x[3], x[4], x[5])
    background_rate = ibd_bkg_rate + accidental_rate
    dwell_time = background_rate/pow(signal_rate/3, 2)
    return dwell_time

bounds = [(0.5,3.5), (10, 60), (0.5, 3.5), (10, 60), (0, 800), (0, 200)]
results = dict()
results['shgo'] = optimize.shgo(dwell_time_opt, bounds)
print(results['shgo'])
'''

print('dwell time (scaling) = ',dwell_time_scaling(singles, 1.2, 8, 1.2, 16, 100, 200, 3, True, True, True),' days\n')
print('dwell time (spectra) = ',dwell_time_heysham(singles, 1.2, 8, 1.2, 16, 100, 200, 3, True, True, True),' days\n')
