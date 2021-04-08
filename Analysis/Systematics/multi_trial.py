import sys
sys.path.append('../')
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




##COMPUTE PRESELECTION-EFFICIENCY
df_rec = uproot.open("../../Tables/SignalTable_20l2_mtexp.root")["SignalTable"].arrays(library="pd")
df_sim = uproot.open("../../Tables/SignalTable_20l2_mtexp.root")["GenTable"].arrays(library="pd")

presel_eff = len(df_rec)/len(df_sim.query("abs(rapidity)<0.5"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
df_mc = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")
score_cuts_array = np.load("../../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
n_events040 = np.sum(uproot.open('../../Utils/AnalysisResults_pPb.root')['AliAnalysisTaskHyperTriton2He3piML_default_summary;1'][11].counts()[:40-1])
branching_ratio = 0.25

signal_list040 = []
error_list040 = []


bkg_models = ['pol0', 'pol1', 'expo']
ff = ROOT.TFile("../../Results/syst_fits.root", "recreate")
bkg_dir = ff.mkdir('bkg_pdf')

for eff,cut in zip(bdt_eff_array, score_cuts_array):
    for bkg_model in bkg_models:
        cut_string = f"model_output>{cut}"
        data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40 and abs(fZ) < 10")["m"])
        mc_data = np.array(df_mc.query(cut_string + " and 2.96<m<3.04")["m"])
        mc_data = mc_data[0:1000]
        res_template = hp.unbinned_mass_fit_mc(data040, eff, 'pol1', mc_data, ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 36, ws_name = f'ws_eff_{eff}')
        signal_list040.append(res_template[0])
        error_list040.append(res_template[1])

ff.cd()

signal_array040 = np.array(signal_list040)
error_array040 = np.array(error_list040)



bdt_eff_array = np.repeat(bdt_eff_array, 3)

    
corrected_counts = signal_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
corrected_error = error_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
Yield = corrected_counts / 2 / 0.98
Yield_error = corrected_error / 2 / 0.98

print('Systemtic uncertainty: ', 100*(np.std(Yield) / 6.27e-07), '%')


numpy_hist = np.histogram(Yield, bins=15)
ghost_hist = aghast.from_numpy(numpy_hist)
root_hist = aghast.to_root(ghost_hist, 'multi_trial')

Fcn1 = ROOT.TLine(6.27e-07, root_hist.GetMinimum(), 6.27e-07, root_hist.GetMaximum())
Fcn1.SetLineWidth(1)
Fcn1.SetLineStyle(2)
Fcn1.SetLineColor(ROOT.kRed - 3)




cv = ROOT.TCanvas(f"syst_multi_trial")

root_hist.Draw()
Fcn1.Draw('lsame')
leg1 = ROOT.TLegend(.68,0.26,.78,0.43)
leg1.SetFillStyle(0)
leg1.SetMargin(0.2)
leg1.SetBorderSize(0)
leg1.SetTextFont(42)
leg1.SetTextSize(0.025)
leg1.AddEntry(root_hist, "Yields from multi-trial", "l")
leg1.AddEntry(Fcn1, "Working point yield", "fl")
leg1.Draw('same')
cv.Write()


root_hist.Write()




histo_barlow = ROOT.TH1D('barlow_check', 'barlow_check', 20, -3 , 5)
wp_index = (bdt_eff_array == 0.72)
wp_yield = Yield[wp_index][0]
wp_error = Yield_error[wp_index][0]


for corr_count, corr_err in zip(Yield,Yield_error):
    if(wp_error!=corr_err):
        correl = (wp_yield - corr_count)/np.sqrt(np.abs(wp_error**2 - corr_err**2))
        histo_barlow.Fill(correl)
histo_barlow.Write()