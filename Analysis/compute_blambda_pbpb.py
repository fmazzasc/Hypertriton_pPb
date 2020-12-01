import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import helpers as hp
import mplhep
import ROOT
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
matplotlib.use("pdf")
matplotlib.style.use(mplhep.style.ALICE)


ml_model_path = "../Utils/Models"

signalH = TreeHandler('../Tables/PbPb/test_set.parquet').apply_preselections('2<pt<8', inplace=False)
model_hdl = ModelHandler()
model_hdl.load_model_handler(ml_model_path + "/model_hndl.pkl")
signalH.apply_model_handler(model_hdl)
df_rec = signalH.get_data_frame()
print("-------------------------------------")


## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
score_cuts_array = np.load("../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../Utils/Efficiencies/bdt_eff_syst_arr.npy")

mean_list = []
mean_error_list = []
mc_corr_list = []


ff = ROOT.TFile("../Results/inv_mass_fits_blam.root", "recreate")
ff.cd()

for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data = np.array(df.query(cut_string + " and 2.96<m<3.04")["m"])
    mc = np.array(df_rec.query(cut_string)["m"])
    mean = np.mean(mc)
    mc_corr_list.append(np.mean(mc))
    plt.figure()
    plt.hist(mc, bins=300)
    plt.plot([mean,mean], [0,30000], 'r--')
    plt.savefig(f'mc_{eff}.png')
    res = hp.unbinned_mass_fit(data, eff, 'pol2', ff, [0,90], [2,8], [0,35], split="", bins = 35)
    mean_list.append(res[4])
    mean_error_list.append(res[5])

ff.Close()
    
mean_array = 1000*np.array(mean_list)
print(mean_array)
print(np.array(mc_corr_list)*1000)
mean_error_array = 1000*np.array(mean_error_list)
mc_corr_array = np.array(mc_corr_list)*1000 - 2991.31
print(mc_corr_array)

b_lam_array = 1115.683 + 1875.61294257 - (mean_array - mc_corr_array)

print("-------------------------------------")
for eff, blam in zip(bdt_eff_array, b_lam_array):
    print('BDT eff: ', eff, '   B_lam: ', blam)
