import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import helpers as hp
import mplhep
import aghast
import ROOT
ROOT.gROOT.SetBatch()

matplotlib.use("pdf")
matplotlib.style.use(mplhep.style.ALICE)


def normalize_ls(data_counts, ls_counts, bins):
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    side_region = np.logical_or(bin_centers<2.992-4*0.0025, bin_centers>2.992+4*0.0025)
    
    side_data_counts = np.sum(data_counts[side_region])
    side_ls_counts = np.sum(ls_counts[side_region])
    scaling_factor = side_data_counts/side_ls_counts
    return scaling_factor

def h1_invmass(counts, bins=34, name=''):
    th1 = ROOT.TH1D(f'{name}', f'{name}_x', int(bins), 2.96, 3.04)

    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
        # th1.SetBinError(index + 1, np.sqrt(counts[index]))

    th1.SetDirectory(0)

    return th1


##COMPUTE PRESELECTION-EFFICIENCY
df_rec = uproot.open("../Tables/SignalTable_pp13TeV_mtexp.root")["SignalTable"].pandas.df().query("pt>0 and rej_accept>0")
df_sim = uproot.open("../Tables/SignalTable_pp13TeV_mtexp.root")["SignalTable"].pandas.df()

presel_eff = len(df_rec)/len(df_sim.query("abs(gY)<0.5 and rej_accept>0"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_data_pp.parquet.gzip")
df_ls = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_ls_pp.parquet.gzip")
df_em = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_em_pp.parquet.gzip")

score_cuts_array = np.load("../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../Utils/Efficiencies/bdt_eff_syst_arr.npy")



ff = ROOT.TFile("../Results/inv_mass_fits_pp_bkg_sub.root", "recreate")
ff.cd()

n_bins = 30

for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data = np.array(df.query(cut_string + " and 2.96<m<3.04")["m"])
    ls_data = np.array(df_ls.query(cut_string + " and 2.96<m<3.04")["m"])
    em_data = np.array(df_em.query(cut_string + " and 2.96<m<3.04")["m"])

    counts_hist_data, bins = np.histogram(data, bins=n_bins)
    counts_hist_ls, _ = np.histogram(ls_data, bins=n_bins)
    counts_hist_ls = counts_hist_ls*normalize_ls(counts_hist_data, counts_hist_ls, bins)

    counts_hist_em, _ = np.histogram(em_data, bins=n_bins)
    counts_hist_em = counts_hist_em*normalize_ls(counts_hist_data, counts_hist_em, bins)

    root_hist_data = h1_invmass(counts_hist_data, n_bins, f"histo_{eff}")
    root_hist_ls = h1_invmass(counts_hist_ls, n_bins, f"histo_{eff}_ls")
    root_hist_em = h1_invmass(counts_hist_em, n_bins, f"histo_{eff}_em")

    root_hist_subtr = h1_invmass(counts_hist_data - counts_hist_ls, n_bins, f"histo_{eff}")
    canv = ROOT.TCanvas(f"data_and_ls_{eff}")
    root_hist_data.Draw("same")
    root_hist_ls.SetMarkerColor(ROOT.kRed)
    root_hist_ls.SetLineColor(ROOT.kRed)
    root_hist_ls.Draw("same")

    # root_hist_em.SetMarkerColor(ROOT.kGreen)
    # root_hist_em.SetLineColor(ROOT.kGreen)
    # root_hist_em.Draw("same")

    canv.Write()

    bkg_params = hp.fit_hist_ls(root_hist_ls, model="pol2", mode=3)

    res_template = hp.fit_hist(root_hist_data, [0,100], [0,9], [0,35], nsigma=3, model=bkg_params, fixsigma=-1, sigma_limits=None, mode=3, split='')

