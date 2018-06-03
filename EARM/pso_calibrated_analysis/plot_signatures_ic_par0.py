from tropical.clustering import PlotSequences, ClusterSequences
import pickle

with open('earm_signatures_ic_par1.pickle', 'rb') as handle:
    all_signatures = pickle.load(handle)

signatures = all_signatures[37]['consumption']
clus = ClusterSequences(seqdata=signatures, unique_sequences=False, truncate_seq=50)

signs_plot = PlotSequences(clus, no_clustering=True)
signs_plot.modal_plot()
