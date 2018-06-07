import numpy as np
from pysb.simulator.scipyode import ScipyOdeSimulator
from earm2_flat import model
import pickle
from pathos.multiprocessing import ProcessingPool as Pool

clus_sp37_labels = np.load('clus_labels_sampled_kd_spectral.npy')

# Use mBid clustering labels
unique_labels = set(clus_sp37_labels)
pars = np.load('sampled_kd_pars.npy')
# pars_ref = pars[0]
tspan = np.linspace(0, 20000, 100)


def sims_kd(label):
    pars_ref1 = np.copy(pars[0])
    pars_ref2 = np.copy(pars[0])

    if label == 0:
        """
        Cluster 0:
        Dominant reactions:
        BidM_BaxC, BidM_Bcl2M, Bid_BclxLM
        """

        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of bax
        pars_label1[:, 58] = pars_ref1[58] * 0.2  # 20% Knock down of bcl2
        # pars_label1[:, 56] = pars_ref[56] * 0.8  # 20% knocl down of bclxl
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8  # 20% Knock down of bax
        pars_label2[:, 57] = pars_ref2[57] * 0.2  # 20% Knock down of mcl1
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 1:
        """
        Cluster 1:
        Dominant reactions:
        BidM_BaxC, BidM_BclxLM, Bid_Mcl1M
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of bax
        # pars_label1[:, 64] = pars_label1[64] * 0.8  # 20% Knock down of bak
        pars_label1[:, 56] = pars_ref1[56] * 0.2  # 20% Knock down of bclxl
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8  # 20% Knock down of bax
        # pars_label2[:, 64] = pars_label2[64] * 0.8  # 20% Knock down of bak
        pars_label2[:, 58] = pars_ref2[58] * 0.2  # 20% Knock down of bcl2
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 2:
        """
        Cluster 2:
        Dominant reactions:
        BidM_BaxC, BidM_Bcl2M, Bid_Mcl1 M
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of bax
        pars_label1[:, 57] = pars_ref1[57] * 0.2  # 20% Knock down of mcl1
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8  # 20% Knock down of bax
        pars_label2[:, 56] = pars_ref2[56] * 0.2  # 20% Knock down of bclxl
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 3:
        """
        Cluster 3:
        Dominant reactions:
        BidM_BaxC, BidM_Mcl1M
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of bax
        pars_label1[:, 57] = pars_ref1[57] * 0.2  # 20% Knock down of mcl1
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8  # 20% Knock down of bax
        pars_label2[:, 58] = pars_ref2[58] * 0.2  # 20% Knock down of bcl2
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 4:
        """
        Cluster 4:
        Dominant reactions:
        BidM_BaxC
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        pars_label1[:, 63] = pars_ref1[63] * 1.2  # 20% overexpression of bax
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        pars_label2[:, 63] = pars_ref2[63] * 0.2  # 20% Knock down of bax
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 5:
        """
        Cluster 5:
        Dominant reactions:
        BidM_BaxC, BidM_BclxLM
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of Bax
        pars_label1[:, 56] = pars_ref1[56] * 0.2  # 20% Knock down of BclxL
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8  # 20% Knock down of Bax
        pars_label2[:, 58] = pars_ref2[58] * 0.2  # 20% Knock down of bcl2
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    if label == 6:
        """
        Cluster 6:
        Dominant reactions:
        BidM_BaxC, BidM_Bcl2M
        """
        pars_label1 = pars[np.where(clus_sp37_labels == label)]
        # pars_label1[:, 63] = pars_ref1[63] * 0.8  # 20% Knock down of bax
        pars_label1[:, 58] = pars_ref1[58] * 0.2  # 20% Knock down of Bcl2

        # pars_label1[:, 64] = pars_label1[64] * 0.8  # 80% Knock down of bak
        # pars_label1[:, 57] = pars_label1[57] * 0.8  # 80% Knock down of mcl1
        sim1 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label1).run()
        sim1.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_good{0}.h5'.format(label))

        pars_label2 = pars[np.where(clus_sp37_labels == label)]
        # pars_label2[:, 63] = pars_ref2[63] * 0.8 # 80% Knock down of bax
        # pars_label2[:, 64] = pars_label2[64] * 0.8  # 80% Knock down of bak
        # pars_label2[:, 58] = pars_label2[58] * 0.8  # 80% Knock down of bcl2
        pars_label2[:, 56] = pars_ref2[56] * 0.2  # 80% Knock down of bclxl
        sim2 = ScipyOdeSimulator(model, tspan=tspan, param_values=pars_label2).run()
        sim2.save('sims_kd80_sensitivities_sampled_kd/earm_scipyode_sims_bad{0}.h5'.format(label))

    return


p = Pool(25)
res = p.amap(sims_kd, unique_labels)
res.get()


