import ROOT
import numpy as np
import matplotlib.pyplot as plt
import random
from math import sqrt
from scipy import optimize

# Script for plotting 2D positron reconstruction efficiency plots

# IBD interaction rates [per second]
rates = {"IBDPositron_LIQUID_ibd_p" : 9.663984059177409e-05, 
         "IBDPositronHeyshamSig_LIQUID_ibd_p_hs" : 9.02e-6,
         "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb" : 2.484e-5} 

p_key = "IBDPositron_LIQUID_ibd_p"
hs_key = "IBDPositronHeyshamSig_LIQUID_ibd_p_hs"
hb_key = "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb"
keys = [p_key, hs_key, hb_key]

# Store types of interactions
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
                       "closest_pmt" : "Closest PMT [m]",
                       "dwall" : "Distance to wall [m]",
                       "reco_dwall" : "Distance to wall [m]"}
        self.values = {}
        self.values['reco_x'] = []
        self.values['reco_y'] = []
        self.values['reco_z'] = []
        self.values['mc_x'] = []
        self.values['mc_y'] = []
        self.values['mc_z'] = []
        self.values['mc_energy'] = []
        self.values['inner_hits'] = []
        self.values['veto_hits'] = []
        self.values['reco_dwall'] = []
        self.values['dwall'] = []
        self.values['n9'] = []
        self.values['good_pos'] = []
        self.values['closest_pmt'] = []
        self.nevents = evts

    # Add a new interaction to the lists
    def add(self, reco_x, reco_y, reco_z, mc_x, mc_y, mc_z, mc_energy,
            inner_hits, veto_hits, n9, good_pos, closest_pmt):
        self.values['reco_x'].append(reco_x)
        self.values['reco_y'].append(reco_y)
        self.values['reco_z'].append(reco_z)
        self.values['mc_x'].append(mc_x)
        self.values['mc_y'].append(mc_y)
        self.values['mc_z'].append(mc_z)
        self.values['mc_energy'].append(mc_energy)
        self.values['inner_hits'].append(inner_hits)
        self.values['veto_hits'].append(veto_hits)
        self.values['n9'].append(n9)
        self.values['good_pos'].append(good_pos)
        self.values['closest_pmt'].append(closest_pmt)
        # True distance to the PMT support
        mc_r = sqrt(mc_x**2 + mc_y**2)
        dwall = 6.7 - max(mc_r, abs(mc_z))
        self.values['dwall'].append(dwall)
        # Reconstructed distance to the PMT support
        reco_r = sqrt(reco_x**2 + reco_y**2)
        reco_dwall = 6.7 - max(reco_r, abs(reco_z))
        self.values['reco_dwall'].append(reco_dwall)

    # Replace all lists with numpy arrays
    def numpify(self):
        for key in self.values:
            self.values[key] = np.array(self.values[key])

    # Calculate the reconstruction efficiency with perfect (truth) reconstruction
    def efficiency(self, wall_cut, energy_cut):
        eff = np.count_nonzero((self.values['dwall'] > wall_cut)
                               & (self.values['mc_energy'] > energy_cut)
                               & (self.values['veto_hits'] < 4))
        eff = eff/self.nevents
        return eff

    # Calculate the reconstruction efficiency with actual reconstruction
    def reco_efficiency(self, wall_cut, energy_cut):
        eff = np.count_nonzero((self.values['closest_pmt'] > wall_cut)
                               & (self.values['n9'] > energy_cut)
                               & (self.values['inner_hits'] > 4)
                               & (self.values['veto_hits'] < 4))
        eff = eff/self.nevents
        return eff


# Make 2D plot of reconstruction efficiency map
def plot_2d(singles, wall_key, energy_key, minmax, bins):
    wall_range = np.linspace(minmax[0][0],minmax[0][1], bins)
    energy_range = np.linspace(minmax[1][0],minmax[1][1], bins)
    
    for key in singles:
        plt.cla()
        plt.clf()
        plot_wall = []
        plot_energy = []
        weights = []
        # Loop over all combinations of dwall and energy cut
        for wall in wall_range:
            for energy in energy_range:
                plot_wall.append(wall)
                plot_energy.append(energy)
                # Calculate either reconstructed or true efficiency
                if wall_key == 'closest_pmt':
                    weights.append(singles[key].reco_efficiency(wall, energy))
                else:
                    #print(wall,energy,singles[key].efficiency(wall, energy))
                    weights.append(singles[key].efficiency(wall, energy))
        # Make a 2D histogram
        plt.hist2d(plot_wall,
                   plot_energy,
                   weights=weights,
                   range=minmax,
                   vmin=0,
                   vmax=0.4,
                   bins=[bins,bins])
        plt.xlabel(singles[key].labels[wall_key])
        plt.ylabel(singles[key].labels[energy_key])
        plt.colorbar()
        plt.savefig('plots/'+key+'_'+wall_key+'_'+energy_key+'_eff')


root_file = ROOT.TFile.Open("gd_water_reco.root")

# Read in events from the tree
singles = {}
for event in root_file.tree:
    key = str(event.int_name)
    if key not in keys:
        continue
    if key not in singles:
        print(key)
        singles[key] = SingleEvent(key, rates[key], event.nevents)
    singles[key].add(event.reco_x/1000., event.reco_y/1000., event.reco_z/1000.,
                     event.mc_x/1000., event.mc_y/1000., event.mc_z/1000.,
                     event.mc_energy, event.inner_hits, event.veto_hits,
                     event.n9, event.good_pos, event.closest_pmt/1000.)

# Turn all lists into numpy arrays for speed
for key in singles:
    singles[key].numpify()

# Make 2D efficiency plots
plot_2d(singles, 'dwall', 'mc_energy', [[0.5, 3.5], [0, 4]], 40)
plot_2d(singles, 'closest_pmt', 'n9', [[0.5, 3.5], [5, 55]], 40)
