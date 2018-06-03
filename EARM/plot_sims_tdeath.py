from earm2_flat import model
from tropical.cluster_analysis import AnalysisCluster
from tropical import util
import numpy as np
import pickle

sims_cupsoda = 'sims_cupsoda.h5'
sims_cupsoda_single = 'sims_cupsoda_single.h5'
sims_scipyode = 'earm_scipyode_sims.h5'
sims_scipyode_single = 'sims_scipyode_single.h5'
sims_cupsoda_sampled_kd = 'earm_cupsoda_sampled_kd.h5'

with open ('sims_sampled_kd_list', 'rb') as fp:
    all_simulations = pickle.load(fp)

tspan = np.linspace(0, 20000, 100)
param_values = np.load('sampled_kd_pars.npy')
simulations = {'trajectories': all_simulations, 'param_values':param_values, 'tspan': tspan}


a = AnalysisCluster(model, clusters=None, sim_results=simulations)

a.plot_dynamics_cluster_types([39], norm=False)
