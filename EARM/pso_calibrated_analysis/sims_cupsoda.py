from pysb.simulator.cupsoda import CupSodaSimulator
from earm2_flat import model
import numpy as np

tspan = np.linspace(0, 20000, 100)

vol = 1e-19

parameters = np.load('calibrated_6572pars.npy')

integrator_opt = {'rtol': 1e-6, 'atol': 1e-6, 'mxsteps': 20000}
cupsoda_solver = CupSodaSimulator(model, tspan, verbose=False, gpu=0, memory_usage='shared_constant', vol=vol, integrator_options=integrator_opt)
sims = cupsoda_solver.run(param_values=parameters)
sims.save("sims_cupsoda.h5")
