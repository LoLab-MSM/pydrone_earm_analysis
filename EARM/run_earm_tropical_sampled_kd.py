import pickle
from tropical.dynamic_signatures_range import run_tropical_multi
from earm2_flat import model
import numpy as np

#with open ('outfile', 'rb') as fp:
#    all_simulations = pickle.load(fp)

with open ('sims_sampled_kd_ic_list.pickle', 'rb') as fp:
    all_simulations = pickle.load(fp)

tspan = np.linspace(0, 20000, 100)
#param_values = np.load('sampled_kd_pars.npy')
param_values = np.load('sampled_kd_ic_parameter.npy')
simulations = {'trajectories': all_simulations, 'param_values':param_values, 'tspan': tspan}

a = run_tropical_multi(model, simulations=simulations, cpu_cores=31, verbose=True)

with open('earm_signatures_sampled_kd_ic.pickle', 'wb') as handle:
    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
