import pickle
from tropical.dynamic_signatures_range import run_tropical_multi
from earm2_flat import model

a = run_tropical_multi(model, simulations='earm_scipyode_sims.h5', cpu_cores=30, verbose=True)

with open('earm_all_kpars_signatures.pickle', 'wb') as handle:
    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
