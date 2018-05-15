from earm2_flat import model
from tropical.cluster_analysis import AnalysisCluster
from tropical import util

sims_cupsoda = 'sims_cupsoda.h5'
sims_scipyode = 'earm_scipyode_sims.h5'

a = AnalysisCluster(model, clusters=None, sim_results=sims_scipyode)

a.plot_dynamics_cluster_types([39], norm=True)
