from earm2_flat import model
from tropical.cluster_analysis import AnalysisCluster
from tropical import util
import pickle 
import numpy as np

sim_0_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good0.h5'
#sim_0_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad0.h5'
sim_1_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good1.h5'
#sim_1_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad1.h5'
sim_2_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good2.h5'
#sim_2_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad2.h5'
sim_3_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good3.h5'
#sim_3_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad3.h5'
sim_4_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good4.h5'
#sim_4_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad4.h5'
sim_5_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good5.h5'
#sim_5_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad5.h5'
sim_6_good = 'sims_baxkd80_sensitivities_sampled_kd/earm_scipyode_sims_good6.h5'
#sim_6_bad = 'sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad6.h5'


with open ('sims_sampled_kd_list', 'rb') as fp:
    all_simulations = pickle.load(fp)

tspan = np.linspace(0, 20000, 100)
#param_values = np.load('sampled_kd_pars.npy')
param_values = np.load('sampled_kd_pars.npy')
simulations_wt = {'trajectories': all_simulations, 'param_values':param_values, 'tspan': tspan}

a = AnalysisCluster(model, clusters=None, sim_results=simulations_wt)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='cluswt', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_0_good)

# a.hist_plot_clusters([82, 83, 84, 85, 86, 87], save_path='figures/')
# a.plot_sp_ic_overlap([82, 83, 84, 85, 86, 87], save_path='figures/')
# a.violin_plot_sps([82, 83, 84, 85, 86, 87], save_path='figures/')
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus0', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_1_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus1', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})


a = AnalysisCluster(model, clusters=None, sim_results=sim_2_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus2', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_3_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus3', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_4_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus4', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_5_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus5', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

a = AnalysisCluster(model, clusters=None, sim_results=sim_6_good)
a.plot_dynamics_cluster_types([39], save_path='sampled_kd_figures/baxkd_80_sensitivities/', fig_label='clus6', species_ftn_fit={39: util.sig_apop},
                              norm=True, **{'p0': [100, 100, 100]})

# a.scatter_plot_pars(ic_par_idxs=[82, 85, 86, 87], cluster=3, save_path='figures/')
