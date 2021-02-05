import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import helpers as hp
import mplhep
import ROOT
ROOT.gROOT.SetBatch()

matplotlib.use("pdf")
matplotlib.style.use(mplhep.style.ALICE)




##COMPUTE PRESELECTION-EFFICIENCY
df_rec = uproot.open("../Tables/SignalTable_pp13TeV_mtexp.root")["SignalTable"].pandas.df().query("pt>0 and bw_accept>0")
df_sim = uproot.open("../Tables/SignalTable_pp13TeV_mtexp.root")["SignalTable"].pandas.df()

presel_eff = len(df_rec)/len(df_sim.query("abs(gY)<0.5 and bw_accept>0"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
df_mc = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_mc_pp.parquet.gzip")
score_cuts_array = np.load("../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
branching_ratio = 0.25

signal_list040 = []
error_list040 = []
delta_mass_list = []
delta_mass_error_list = []
mean_mass_list = []

ff = ROOT.TFile("../Results/inv_mass_fits_pp.root", "recreate")
bkg_dir = ff.mkdir('bkg_pdf')
ff.cd()


for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40")["m"])
    mc_data = np.array(df_mc.query(cut_string + " and 2.96<m<3.04")["m"])
    mean_mass_list.append(np.mean(mc_data))
    mc_data = mc_data[0:20000]
    res_template = hp.unbinned_mass_fit(data040, eff, 'pol1', ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 34)

    signal_list040.append(res_template[0])
    error_list040.append(res_template[1])
    delta_mass_list.append(res_template[-2])
    delta_mass_error_list.append(res_template[-1])


signal_array040 = np.array(signal_list040)
error_array040 = np.array(error_list040)
mean_mass_array040 = np.array(mean_mass_list)
delta_mass_array = np.array(delta_mass_list)
delta_mass_error_array = np.array(delta_mass_error_list)