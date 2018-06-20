import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pathos.multiprocessing import ProcessingPool as Pool
from earm2_flat import model
import pickle

tspan = np.linspace(0, 20000, 100)


def run_simulation(param_values):
    sim = ScipyOdeSimulator(model, tspan=tspan).run(param_values=param_values).species
    return sim

all_parameters = np.load('sampled_kd_ic_parameter.npy')
cpu_cores = 31
p = Pool(cpu_cores)
res = p.amap(run_simulation, all_parameters)
sims = res.get()

with open('sims_sampled_kd_ic_list.pickle', 'wb') as fp:
    pickle.dump(sims, fp)
