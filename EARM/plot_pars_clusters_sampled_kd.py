from earm2_flat import model
from tropical.cluster_analysis import AnalysisCluster
from tropical import util
import numpy as np
import pickle

with open ('sims_sampled_kd_list', 'rb') as fp:
    all_simulations = pickle.load(fp)

tspan = np.linspace(0, 20000, 100)
param_values = np.load('sampled_kd_pars.npy')
simulations = {'trajectories': all_simulations, 'param_values':param_values, 'tspan': tspan}

cluster_labels = np.load('clus_labels_sampled_kd.npy')
a = AnalysisCluster(model, clusters=cluster_labels, sim_results=simulations)

a.plot_dynamics_cluster_types([37], norm=False, save_path='sampled_kd_figures')
a.violin_plot_sps([82, 86, 84, 70, 73, 96, 98, 100, 88, 92, 94], save_path='sampled_kd_figures')
a.violin_plot_kd([(83, 82), (87, 86), (85, 84), (71, 70), (74, 73), (97, 96), (99, 98), (101, 100), (89, 88), (91, 90), (93, 92), (95, 94)], save_path='sampled_kd_figures')
