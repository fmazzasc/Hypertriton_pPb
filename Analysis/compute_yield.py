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
df_rec = uproot.open("../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df()
df_sim = uproot.open("../Tables/SignalTable_17d_mtexp.root")["GenTable"].pandas.df()


presel_eff = len(df_rec)/len(df_sim.query("abs(rapidity)<0.5"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
df_mc = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")
score_cuts_array = np.load("../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
n_events040 = np.sum(uproot.open('../Utils/AnalysisResults_pPb.root')['AliAnalysisTaskHyperTriton2He3piML_default_summary;1'][11].values[:40-1])
branching_ratio = 0.25

signal_list040 = []
error_list040 = []
n1_list040 = []


ff = ROOT.TFile("../Results/inv_mass_fits.root", "recreate")
bkg_dir = ff.mkdir('bkg_pdf')
ff.cd()


for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40 and abs(fZ)<10")["m"])
    res040 = hp.unbinned_mass_fit(data040, eff, 'pol0', ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 34)

    mc_data = np.array(df_mc.query(cut_string + " and 2.96<m<3.04")["m"])
    mean_mass = np.mean(mc_data)
    hist = ROOT.TH1D(f'histo_mc_{eff}', f'histo_mc_{eff}', 500, 2.96, 3.04)
    for mc_entry in mc_data:
        hist.Fill(mc_entry)# - mean_mass + res040[4])

    res_template = hp.unbinned_mass_fit_mc(data040, eff, 'pol1', hist, ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 34, sign_range = [res040[-2], res040[-1]], ws_name=f"ws_pol1_{eff*100:.0f}")

    signal_list040.append(res_template[0])
    error_list040.append(res_template[1])
    n1_list040.append(res_template[2])


signal_array040 = np.array(signal_list040)
error_array040 = np.array(error_list040)


np.save('../Utils/pdf_fraction_040', n1_list040)



###COMPUTE YIELD############################


corrected_counts = signal_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
print(corrected_counts/2)
corrected_error = error_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
Yield = np.float64(corrected_counts[bdt_eff_array==selected_bdt_eff]  / 2)
Yield = Yield/0.98   #absorption correction

stat_error = np.float64(corrected_error[bdt_eff_array==selected_bdt_eff] / 2)
syst_error = float(np.std(corrected_counts) / corrected_counts[bdt_eff_array==selected_bdt_eff])
syst_error = syst_error*Yield
pt_shape_syst = 0.00628*Yield
abs_syst = 0.041*Yield
fit_syst = 0.055*Yield
syst_error = np.sqrt(syst_error**2 + pt_shape_syst**2 + abs_syst**2 + fit_syst**2)

print("-------------------------------------")
print(f"Yield [(matter + antimatter) / 2] 0-40% = {Yield:.3e} +- {stat_error:.3e} (stat.) +- {syst_error:.3e} (syst.)")
#############################################

###COMPUTE S3
protonVals = [{"bin":[0.0, 5.0], "measure": [2.280446e+00, 4.394273e-03, 1.585968e-01, 7.001060e-02]},
{"bin":[5.0, 10.0], "measure": [1.853576e+00, 3.901768e-03, 1.284426e-01, 5.548104e-02]},
{"bin":[10.0, 20.0], "measure": [1.577375e+00, 2.739986e-03, 1.084793e-01, 4.467687e-02]},
{"bin":[20.0, 40.0], "measure": [1.221875e+00, 1.735134e-03, 8.190972e-02, 3.347605e-02]},
{"bin":[40.0, 60.0], "measure": [8.621290e-01, 1.381181e-03, 5.608663e-02, 2.265973e-02]},
{"bin":[60.0, 80.0], "measure": [5.341445e-01, 1.034190e-03, 3.408551e-02, 1.377331e-02]},
{"bin":[80.0, 100.0], "measure": [2.307767e-01, 6.969543e-04, 1.450569e-02, 6.911692e-03]}]
lambdaVals = [{"bin": [0.0, 5.0], "measure": [1.630610e+00, 7.069538e-03, 1.448270e-01, 7.144051e-02]},
{"bin": [5.0, 10.0], "measure": [1.303508e+00, 6.313141e-03, 1.122199e-01, 5.217514e-02]},
{"bin": [10.0, 20.0], "measure": [1.093400e+00, 4.432728e-03, 9.410049e-02, 4.447422e-02]},
{"bin": [20.0, 40.0], "measure": [8.332270e-01, 3.036428e-03, 6.919569e-02, 3.205334e-02]},
{"bin": [40.0, 60.0], "measure": [5.678968e-01, 2.357040e-03, 4.722551e-02, 2.110751e-02]},
{"bin": [60.0, 80.0], "measure": [3.355044e-01, 1.855427e-03, 2.871323e-02, 1.389226e-02]},
{"bin": [80.0, 100.0], "measure": [1.335216e-01, 1.353847e-03, 1.142858e-02, 6.444149e-03]}]



protonAv040 = hp.computeAverage(protonVals,40)
lambdaAv040 = hp.computeAverage(lambdaVals,40)
He3Av040 = [2.07e-6, 1.05e-7, 1.95e-7]


# d$N$/d$\eta$ obtained by simple weighted average of the values published in https://arxiv.org/pdf/1910.14401.pdf

x_pPb040=np.array([29.4], dtype=np.float64)
xe_pPb040=np.array([0.6], dtype=np.float64)


s3_040 = np.array([Yield * protonAv040[0] / (He3Av040[0] * lambdaAv040[0])], dtype=np.float64)
s3stat040 = np.array([s3_040[0] * hp.myHypot(stat_error / Yield, He3Av040[1] / He3Av040[0], protonAv040[1] / protonAv040[0], lambdaAv040[1] / lambdaAv040[0])], dtype=np.float64)
s3syst040 = np.array([s3_040[0] * hp.myHypot(syst_error / Yield, He3Av040[2] / He3Av040[0], protonAv040[2] / protonAv040[0], lambdaAv040[2] / lambdaAv040[0])], dtype=np.float64)

print("-------------------------------------")
print(f"S3 0-40% = {s3_040[0]:.2e} +- {s3stat040[0]:.2e} (stat.) +- {s3syst040[0]:.2e} (syst.)")

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

s3_csm = ROOT.TGraphErrors("../Utils/ProdModels/FullCSM-S3.dat","%lg %*s %*s %*s %*s %*s %*s %*s %*s %*s %lg")
fin = ROOT.TFile('../Utils/ProdModels/CSM_predictions_S3_T155MeV_Vc3dNdy.root')
s3_csm_std = fin.Get('gCSM_S3_vs_dNchdEta')


s3_csm.SetLineColor(kOrangeC)
s3_csm.SetLineWidth(1)
s3_csm.SetTitle("Full canonical SHM")

s3_csm_std.SetLineColor(kGreenC)
s3_csm_std.SetLineWidth(1)
s3_csm_std.SetTitle("Canonical SHM with T=155MeV, Vc = 3dV/dy")


s3_2body = ROOT.TGraphErrors("../Utils/ProdModels/s3_2body.csv","%lg %lg %*s","\t,")
s3_2body.SetLineColor(kBlueC)
s3_2body.SetMarkerColor(kBlueC)
s3_2body.SetTitle("2-body coalescence")
s3_3body = ROOT.TGraphErrors("../Utils/ProdModels/s3_3body.csv","%lg %lg %*s","\t,")
s3_3body.SetLineColor(kAzureC)
s3_3body.SetMarkerColor(kAzureC)
s3_3body.SetTitle("3-body coalescence")

cv = ROOT.TCanvas("cv")
cv.SetBottomMargin(0.14)
frame=cv.DrawFrame(9,0.01,2200,1.0,";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};S_{3}")
frame.GetXaxis().SetTitleOffset(1.25)
cv.SetLogx()
# cv.SetLogy()
s3_csm.Draw("L")
s3_csm_std.Draw('L')
s3_2body.Draw("L")
s3_3body.Draw("L")

x = np.array([1447], dtype=np.float64)
ex = np.array([39], dtype=np.float64)
y = np.array([0.6], dtype=np.float64)
ey = np.array([0.13], dtype=np.float64)
eys = np.array([0.21], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_stat.SetLineColor(kGreyC)
pbpb_stat.SetMarkerColor(kGreyC)
pbpb_stat.SetMarkerStyle(20)
pbpb_stat.Draw("Pz")


pbpb_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_syst.SetTitle("ALICE Pb-Pb #sqrt{#it{s}_{NN}}=2.76 TeV")
pbpb_syst.SetLineColor(kGreyC)
pbpb_syst.SetMarkerColor(kGreyC)
pbpb_syst.SetFillStyle(0)
pbpb_syst.SetMarkerStyle(20)
pbpb_syst.Draw("P2")




ppb_stat040 = ROOT.TGraphErrors(1,x_pPb040,s3_040,zero,s3stat040)
ppb_stat040.SetLineColor(kVioletC)
ppb_stat040.SetMarkerColor(kVioletC)
ppb_stat040.SetMarkerStyle(20)

ppb_syst040 = ROOT.TGraphErrors(1,x_pPb040, s3_040, xe_pPb040, s3syst040)
ppb_syst040.SetTitle("ALICE Internal p-Pb, 0-40%, #sqrt{#it{s}_{NN}}=5.02 TeV")
ppb_syst040.SetLineColor(kVioletC)
ppb_syst040.SetMarkerColor(kVioletC)
ppb_syst040.SetFillStyle(0)
ppb_syst040.SetMarkerStyle(20)


leg = ROOT.TLegend(0.15,0.75,0.7,0.85)
leg.SetMargin(0.14)


ppb_stat040.Draw("Pz")
ppb_syst040.Draw("P2")
leg.AddEntry(ppb_syst040,"","pf")

pinfo = ROOT.TPaveText(0.18,0.68,0.33,0.73, 'NDC')
pinfo.SetBorderSize(0)
pinfo.SetFillStyle(0)
pinfo.SetTextAlign(30+3)
pinfo.SetTextFont(42)
pinfo.AddText('B.R. = 0.25')
pinfo.Draw()


leg.AddEntry(pbpb_syst,"","pf")
leg.SetEntrySeparation(0.2)
legT = ROOT.TLegend(0.58,0.15,0.92,0.47)
legT.SetEntrySeparation(0.2)
legT.SetMargin(0.14)
legT.AddEntry(s3_csm)
legT.AddEntry(s3_csm_std)
legT.AddEntry(s3_2body)
legT.AddEntry(s3_3body)
leg.SetFillStyle(0)
legT.SetFillStyle(0)
leg.Draw()
legT.Draw()

cv.Draw()

cv.SaveAs("../Results/s3.pdf")
cv.SaveAs("../Results/s3.png")
