import sys
sys.path.append('../')
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


parser = argparse.ArgumentParser()
parser.add_argument("-c", '--cent040', help='Plot only 0-40% point', action = 'store_true')
parser.add_argument("-b", '--both', help='Plot both points', action = 'store_true')
args = parser.parse_args()

only040 = args.cent040
onlyInt = not only040
if args.both:
    onlyO40 = False
    onlyInt = False





##COMPUTE PRESELECTION-EFFICIENCY
df_rec = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df()
df_sim = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["GenTable"].pandas.df()

presel_eff = len(df_rec)/len(df_sim.query("abs(rapidity)<0.5"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
score_cuts_array = np.load("../../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
signal_list = []
error_list = []
signal_list040 = []
error_list040 = []



n1_list040 = []
n1_list = []

ff = ROOT.TFile("../../Results/inv_mass_fits.root", "recreate")
bkg_dir = ff.mkdir('bkg_pdf')
ff.cd()

for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data = np.array(df.query(cut_string + " and 2.96<m<3.04 ")["m"])
    data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40 and abs(fZ) < 10")["m"])
    res = hp.unbinned_mass_fit(data, eff, 'pol0', ff, bkg_dir, [0,90], [0,10], [0,35], split="", bins = 34)
    res040 = hp.unbinned_mass_fit(data040, eff, 'pol0', ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 34)
    signal_list.append(res[0])
    error_list.append(res[1])
    signal_list040.append(res040[0])
    error_list040.append(res040[1])
    n1_list.append(res[-1])
    n1_list040.append(res040[-1])


signal_array = np.array(signal_list)
error_array = np.array(error_list)
signal_array040 = np.array(signal_list040)
error_array040 = np.array(error_list040)

signal_array_list = [signal_array, signal_array040]
error_array_list = [error_array, error_array040]

np.save('../../Utils/pdf_fraction_040', n1_list040)
np.save('../../Utils/pdf_fraction_040', n1_list040)

###COMPUTE YIELD
n_events = 7.38470e+08 + 1.254483e+8
n_events040 = np.sum(uproot.open('../../Utils/AnalysisResults_pPb.root')['AliAnalysisTaskHyperTriton2He3piML_default_summary;1'][11].values[:40-1])
n_events_list = [n_events, n_events040]
string_list = ['', 'in 0-40%']
branching_ratio = 0.25

yield_list = []
stat_list = []
syst_list = []

for sig_arr, err_arr, n_ev, string in zip(signal_array_list, error_array_list, n_events_list, string_list):
    corrected_counts = sig_arr/n_ev/branching_ratio/presel_eff/bdt_eff_array
    corrected_error = err_arr/n_ev/branching_ratio/presel_eff/bdt_eff_array
    Yield = np.float64(corrected_counts[bdt_eff_array==selected_bdt_eff]  / 2)
    stat_error = np.float64(corrected_error[bdt_eff_array==selected_bdt_eff] / 2)
    syst_error = np.std(corrected_counts) / 2
    pt_shape_syst = 0.0628*Yield
    abs_syst = 0.0577*Yield
    syst_error += pt_shape_syst + abs_syst 
    yield_list.append(Yield)
    stat_list.append(stat_error)
    syst_list.append(syst_error)
    print("-------------------------------------")
    print(f"Yield [(matter + antimatter) / 2] {string} = {Yield:.2e} +- {stat_error:.2e} (stat.) +- {syst_error:.2e} (syst.)")


###COMPUTE hp_ratio
protonVals = [{"bin":[0.0, 5.0], "measure": [2.280446e+00, 4.394273e-03, 1.585968e-01, 7.001060e-02]},
{"bin":[5.0, 10.0], "measure": [1.853576e+00, 3.901768e-03, 1.284426e-01, 5.548104e-02]},
{"bin":[10.0, 20.0], "measure": [1.577375e+00, 2.739986e-03, 1.084793e-01, 4.467687e-02]},
{"bin":[20.0, 40.0], "measure": [1.221875e+00, 1.735134e-03, 8.190972e-02, 3.347605e-02]},
{"bin":[40.0, 60.0], "measure": [8.621290e-01, 1.381181e-03, 5.608663e-02, 2.265973e-02]},
{"bin":[60.0, 80.0], "measure": [5.341445e-01, 1.034190e-03, 3.408551e-02, 1.377331e-02]},
{"bin":[80.0, 100.0], "measure": [2.307767e-01, 6.969543e-04, 1.450569e-02, 6.911692e-03]}]


protonAv = hp.computeAverage(protonVals)
protonAv040 = hp.computeAverage(protonVals,40)

print(protonAv)
print(protonAv040)


# d$N$/d$\eta$ obtained by simple weighted average of the values published in https://arxiv.org/pdf/1910.14401.pdf
x_pPb = np.array([17.8], dtype=np.float64)
x_pPb040=np.array([29.4], dtype=np.float64)
xe_pPb = np.array([0.4], dtype=np.float64)
xe_pPb040=np.array([0.6], dtype=np.float64)

hp_ratio = np.array([yield_list[0] / protonAv[0]], dtype=np.float64)
hp_ratiostat = np.array([hp_ratio[0] * hp.myHypot(stat_list[0] / yield_list[0], protonAv[1] / protonAv[0])], dtype=np.float64)
hp_ratiosyst = np.array([hp_ratio[0] * hp.myHypot(syst_list[0] / yield_list[0], protonAv[2] / protonAv[0])], dtype=np.float64)

hp_ratio_040 = np.array([yield_list[1] / protonAv040[0]], dtype=np.float64)
hp_ratiostat040 = np.array([hp_ratio_040[0] * hp.myHypot(stat_list[1] / yield_list[1], protonAv040[1] / protonAv040[0])], dtype=np.float64)
hp_ratiosyst040 = np.array([hp_ratio_040[0] * hp.myHypot(syst_list[1] / yield_list[1], protonAv040[2] / protonAv040[0])], dtype=np.float64)

print("-------------------------------------")
print(f"Hyp/Proton {string_list[0]} = {hp_ratio[0]:.2e} +- {hp_ratiostat[0]:.2e} (stat.) +- {hp_ratiosyst[0]:.2e} (syst.)")
print(f"Hyp/Proton {string_list[1]} = {hp_ratio_040[0]:.2e} +- {hp_ratiostat040[0]:.2e} (stat.) +- {hp_ratiosyst040[0]:.2e} (syst.)")

kBlueC  = ROOT.TColor.GetColor("#2077b4");
kRedC  = ROOT.TColor.GetColor("#d62827");
kGreenC  = ROOT.TColor.GetColor("#2ba02b");
kOrangeC  = ROOT.TColor.GetColor("#ff7f0f");
kVioletC  = ROOT.TColor.GetColor("#9467bd");
kPinkC  = ROOT.TColor.GetColor("#e377c1");
kGreyC  = ROOT.TColor.GetColor("#7f7f7f");
kBrownC  = ROOT.TColor.GetColor("#8c564c");
kAzureC  = ROOT.TColor.GetColor("#18becf");
kGreenBC  = ROOT.TColor.GetColor("#bcbd21");

hp_ratio_csm = ROOT.TGraphErrors("../../Utils/ProdModels/hyp_p_ratio.dat","%lg %*s %*s %*s %*s %lg %*s")
# fin = ROOT.TFile('../../../Desktop/vanilla_CSM_predictions_H3L_to_P.root')
# hp_ratio_csm = fin.Get('gCSM_3HL_over_p')
print(hp_ratio_csm)
hp_ratio_csm.Draw()
hp_ratio_csm.SetLineColor(kOrangeC)
hp_ratio_csm.SetLineWidth(1)
hp_ratio_csm.SetTitle("Full canonical SHM")
# hp_ratio_2body = ROOT.TGraphErrors("../../Utils/ProdModels/hp_ratio_2body.csv","%lg %lg %*s","\t,")
# hp_ratio_2body.SetLineColor(kBlueC)
# hp_ratio_2body.SetMarkerColor(kBlueC)
# hp_ratio_2body.SetTitle("2-body coalescence")
# hp_ratio_3body = ROOT.TGraphErrors("../../Utils/ProdModels/hp_ratio_3body.csv","%lg %lg %*s","\t,")
# hp_ratio_3body.SetLineColor(kAzureC)
# hp_ratio_3body.SetMarkerColor(kAzureC)
# hp_ratio_3body.SetTitle("3-body coalescence")

cv = ROOT.TCanvas("cv")
cv.SetBottomMargin(0.14)
frame=cv.DrawFrame(5., 1e-7, 200, hp_ratio[0]*34,";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}; _{#Lambda}^{3}H/p")
frame.GetXaxis().SetTitleOffset(1.25)
cv.SetLogx()
cv.SetLogy()
hp_ratio_csm.Draw("L")
# hp_ratio_2body.Draw("L")
# hp_ratio_3body.Draw("L")

zero = np.array([0], dtype=np.float64)




ppb_stat = ROOT.TGraphErrors(1,x_pPb,hp_ratio,zero,hp_ratiostat)
ppb_stat.SetLineColor(kRedC)
ppb_stat.SetMarkerColor(kRedC)
ppb_stat.SetMarkerStyle(20)

ppb_syst = ROOT.TGraphErrors(1,x_pPb,hp_ratio,xe_pPb,hp_ratiosyst)
ppb_syst.SetTitle("ALICE Internal p-Pb, #sqrt{#it{s}_{NN}}=5.02 TeV")
ppb_syst.SetLineColor(kRedC)
ppb_syst.SetMarkerColor(kRedC)
ppb_syst.SetFillStyle(0)
ppb_syst.SetMarkerStyle(20)


ppb_stat040 = ROOT.TGraphErrors(1,x_pPb040,hp_ratio_040,zero,hp_ratiostat040)
ppb_stat040.SetLineColor(kVioletC)
ppb_stat040.SetMarkerColor(kVioletC)
ppb_stat040.SetMarkerStyle(20)

ppb_syst040 = ROOT.TGraphErrors(1,x_pPb040, hp_ratio_040, xe_pPb040, hp_ratiosyst040)
ppb_syst040.SetTitle("ALICE Internal p-Pb, 0-40%, #sqrt{#it{s}_{NN}}=5.02 TeV")
ppb_syst040.SetLineColor(kVioletC)
ppb_syst040.SetMarkerColor(kVioletC)
ppb_syst040.SetFillStyle(0)
ppb_syst040.SetMarkerStyle(20)


leg = ROOT.TLegend(0.55,0.4,0.88,0.5)
leg.SetMargin(0.14)

if only040:
    ppb_stat040.Draw("Pz")
    ppb_syst040.Draw("P2")
    leg.AddEntry(ppb_syst040,"","pf")

elif onlyInt:
    ppb_stat.Draw("Pz")
    ppb_syst.Draw("P2")
    leg.AddEntry(ppb_syst,"","pf")

else:
    ppb_stat040.Draw("Pz")
    ppb_syst040.Draw("P2")    
    ppb_stat.Draw("Pz")
    ppb_syst.Draw("P2")
    leg.AddEntry(ppb_syst,"","pf")
    leg.AddEntry(ppb_syst040,"","pf")


leg.SetEntrySeparation(0.2)
legT = ROOT.TLegend(0.55,0.18,0.88,0.38)
legT.SetMargin(0.14)
legT.AddEntry(hp_ratio_csm)
# legT.AddEntry(hp_ratio_2body)
# legT.AddEntry(hp_ratio_3body)
leg.SetFillStyle(0)
legT.SetFillStyle(0)
leg.Draw()
legT.Draw()

cv.Draw()

cv.SaveAs("../../Results/hp_ratio.pdf")
cv.SaveAs("../../Results/hp_ratio.png")
