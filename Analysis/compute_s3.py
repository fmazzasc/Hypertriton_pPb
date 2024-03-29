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
df_rec = uproot.open("../Tables/SignalTable_20l2_mtexp.root")["SignalTable"].arrays(library="pd")
df_sim = uproot.open("../Tables/SignalTable_20l2_mtexp.root")["GenTable"].arrays(library="pd")

presel_eff = len(df_rec)/len(df_sim.query("abs(rapidity)<0.5"))
print("-------------------------------------")
print("Pre-selection efficiency: ", presel_eff)

## FIT INVARIANT MASS SPECTRA
df = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
df_mc = pd.read_parquet("../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")
score_cuts_array = np.load("../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
n_events040 = np.sum(uproot.open('../Utils/AnalysisResults_pPb.root')['AliAnalysisTaskHyperTriton2He3piML_default_summary;1'][11].counts()[:40-1])
branching_ratio = 0.25

signal_list040 = []
error_list040 = []
delta_mass_list = []
delta_mass_error_list = []
mean_mass_list = []

ff = ROOT.TFile("../Results/inv_mass_fits.root", "recreate")
bkg_dir = ff.mkdir('bkg_pdf')
ff.cd()


for eff,cut in zip(bdt_eff_array, score_cuts_array):
    cut_string = f"model_output>{cut}"
    data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40 and abs(fZ)<10")["m"])
    mc_data = np.array(df_mc.query(cut_string + " and 2.96<m<3.04")["m"])
    mean_mass_list.append(np.mean(mc_data))
    mc_data = mc_data[0:1000]
    res_template = hp.unbinned_mass_fit_mc(data040, eff, 'pol1', mc_data, ff, bkg_dir, [0,40], [0,10], [0,35], split="", cent_string='040', bins = 34, ws_name = f'ws_eff_{eff}')

    signal_list040.append(res_template[0])
    error_list040.append(res_template[1])
    delta_mass_list.append(res_template[-2])
    delta_mass_error_list.append(res_template[-1])


signal_array040 = np.array(signal_list040)
error_array040 = np.array(error_list040)
mean_mass_array040 = np.array(mean_mass_list)
delta_mass_array = np.array(delta_mass_list)
delta_mass_error_array = np.array(delta_mass_error_list)
# print(delta_mass_array)
# B_lam = 1.115683 + 1.87561294257 - (2.99131 - 1000*delta_mass_array[bdt_eff_array==selected_bdt_eff])
# B_lam_error = 1000*delta_mass_error_array[bdt_eff_array==selected_bdt_eff]
# print(f"BLam: {B_lam*1000} +- {B_lam_error*1000}" )


###COMPUTE YIELD############################


corrected_counts = signal_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
corrected_error = error_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
Yield = np.float64(corrected_counts[bdt_eff_array==selected_bdt_eff]  / 2)
Yield = Yield/0.97   #absorption correction


stat_error = np.float64(corrected_error[bdt_eff_array==selected_bdt_eff] / 2)

syst_error = 0.14*Yield
pt_shape_syst = 0.07*Yield
abs_syst = 0.041*Yield
br_syst = 0.09*Yield
# fit_syst = 0.055*Yield
syst_error = np.sqrt(syst_error**2 + pt_shape_syst**2 + abs_syst**2 + br_syst**2)

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
{"bin": [40.0, 60.0], "measure": [5.678968e-01, 2.3570100e-03, 4.722551e-02, 2.110751e-02]},
{"bin": [60.0, 80.0], "measure": [3.355044e-01, 1.855427e-03, 2.871323e-02, 1.389226e-02]},
{"bin": [80.0, 100.0], "measure": [1.335216e-01, 1.353847e-03, 1.142858e-02, 6.444149e-03]}]



protonAv040 = hp.computeAverage(protonVals,40)
lambdaAv040 = hp.computeAverage(lambdaVals,40)


He3Av040 = [2.07e-6, 1.05e-7, 1.95e-7]
print("s---------------------------------------")
print('Proton: ', np.array(protonAv040)/2)
print('Lambda: ', np.array(lambdaAv040)/2)
print('He3: ', He3Av040)
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


s3_csm_1 = ROOT.TGraphErrors("../Utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %*s %lg")
s3_csm_1.SetLineColor(922)
s3_csm_1.SetLineWidth(2)
s3_csm_1.SetTitle("SHM, #it{Vc} = d#it{V}/d#it{y}")

s3_csm_3 = ROOT.TGraphErrors("../Utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.3dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %*s %lg")
s3_csm_3.SetLineColor(922)
s3_csm_3.SetLineWidth(2)
s3_csm_3.SetLineStyle(2)
s3_csm_3.SetTitle("SHM, #it{Vc} = 3d#it{V}/d#it{y}")

n = s3_csm_1.GetN()
grshade = ROOT.TGraph(2*n)
for i in range(n) : 
   grshade.SetPoint(i, s3_csm_3.GetPointX(i), s3_csm_3.GetPointY(i))
   grshade.SetPoint(n + i, s3_csm_1.GetPointX(n - i -1), s3_csm_1.GetPointY(n - i - 1))
   
grshade.SetFillColorAlpha(16, 0.571)





s3_2body = ROOT.TGraphErrors("../Utils/ProdModels/coalescence/s3_2body.csv","%lg %lg %lg")
s3_2body.SetLineColor(kBlueC)
s3_2body.SetMarkerColor(kBlueC)
s3_2body.SetFillStyle(3145)
s3_2body.SetTitle("2-body coalescence")
s3_2body.SetLineWidth(2)

s3_3body = ROOT.TGraphErrors("../Utils/ProdModels/coalescence/s3_3body.csv","%lg %lg %lg")
s3_3body.SetLineColor(kAzureC)
s3_3body.SetMarkerColor(kAzureC)
s3_3body.SetTitle("3-body coalescence")
s3_3body.SetFillStyle(3014)
s3_3body.SetLineWidth(2)




s3_3body.SetFillColorAlpha(kAzureC, 0.571)
s3_2body.SetFillColorAlpha(kBlueC, 0.571)

mg = ROOT.TMultiGraph()
mg.Add(s3_2body)
mg.Add(s3_3body)

cv = ROOT.TCanvas("cv", "cv", 700,700)
cv.SetBottomMargin(0.145)
cv.SetLeftMargin(0.14)
cv.SetTopMargin(0.01)
cv.SetRightMargin(0.01)
frame=cv.DrawFrame(9,0, 3e3,1.0,";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5};#it{S_{3}}")
cv.SetLogx()
mg.Draw("4al same")
# grshade.Draw("f same")
s3_csm_1.Draw("L same")
s3_csm_3.Draw("L same")
mg.GetYaxis().SetRangeUser(0.01, 1)
mg.GetXaxis().SetTitle('#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}')
mg.GetYaxis().SetTitle('#it{S}_{3}')
mg.GetXaxis().SetTitleOffset(1.1)
mg.GetYaxis().SetTitleOffset(0.95)
mg.GetYaxis().SetTitleSize(0.06)
mg.GetXaxis().SetTitleSize(0.06)
mg.GetYaxis().SetLabelSize(0.035)
mg.GetXaxis().SetLabelSize(0.04)
mg.GetYaxis().SetLabelOffset(0.015)


mg.GetYaxis().SetRangeUser(0, 0.99)
mg.GetXaxis().SetRangeUser(5, 3e3)





x = np.array([30.81], dtype=np.float64)
ex = np.array([0.44], dtype=np.float64)
y = np.array([0.195], dtype=np.float64)
ey = np.array([0.049], dtype=np.float64)
eys = np.array([0.044], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pp_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pp_stat.SetLineColor(kOrangeC)
pp_stat.SetMarkerColor(kOrangeC)
pp_stat.SetMarkerStyle(21)
pp_stat.SetMarkerSize(1)
pp_stat.SetLineWidth(1)
# pp_stat.Draw("Pz")


pp_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pp_syst.SetTitle("ALICE Preliminary pp, HM trigger, #sqrt{#it{s}} = 13 TeV")
pp_syst.SetLineColor(kOrangeC)
pp_syst.SetMarkerColor(kOrangeC)
pp_syst.SetFillStyle(0)
pp_syst.SetMarkerStyle(21)
pp_syst.SetMarkerSize(1)
pp_syst.SetLineWidth(1)
# pp_syst.Draw("P2")






x = np.array([1447], dtype=np.float64)
ex = np.array([39], dtype=np.float64)
y = np.array([0.6], dtype=np.float64)
ey = np.array([0.13], dtype=np.float64)
eys = np.array([0.21], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_stat.SetLineColor(ROOT.kBlue + 3)
pbpb_stat.SetMarkerColor(ROOT.kBlue + 3)
pbpb_stat.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_stat.SetMarkerSize(1.5)
pbpb_stat.Draw("Pz")


pbpb_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_syst.SetTitle("ALICE Pb#minusPb, 0#minus10%, #sqrt{#it{s}_{NN}} = 2.76 TeV")
pbpb_syst.SetLineColor(ROOT.kBlue + 3)
pbpb_syst.SetMarkerColor(ROOT.kBlue + 3)
pbpb_syst.SetFillStyle(0)
pbpb_syst.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_syst.SetMarkerSize(1.5)

pbpb_syst.Draw("P2")




ppb_stat040 = ROOT.TGraphErrors(1,x_pPb040,s3_040,zero,s3stat040)
ppb_stat040.SetLineColor(kRedC)
ppb_stat040.SetMarkerColor(kRedC)
ppb_stat040.SetMarkerStyle(20)

ppb_syst040 = ROOT.TGraphErrors(1,x_pPb040, s3_040, xe_pPb040, s3syst040)
ppb_syst040.SetTitle("ALICE p#minusPb, 0#minus40%, #sqrt{#it{s}_{NN}} = 5.02 TeV")
ppb_syst040.SetLineColor(kRedC)
ppb_syst040.SetMarkerColor(kRedC)
ppb_syst040.SetFillStyle(0)
ppb_syst040.SetMarkerStyle(20)





pinfo = ROOT.TPaveText(0.4,0.7, 0.73, 0.76, 'NDC')
pinfo.SetBorderSize(0)
pinfo.SetFillStyle(0)
pinfo.SetTextAlign(30+3)
pinfo.SetTextFont(42)
pinfo.AddText('#it{S_{3}} = ({}_{#Lambda}^{3}H / ^{3}He ) / ( #Lambda / p )')
pinfo.Draw()



pinfo2 = ROOT.TPaveText(0.215,0.865, 0.465, 0.9, 'NDC')
pinfo2.SetBorderSize(0)
pinfo2.SetFillStyle(0)
# pinfo2.SetTextAlign(30+3)
pinfo2.SetTextFont(42)
pinfo2.AddText('B.R. = 0.25 #pm 0.02')
pinfo2.Draw()


leg = ROOT.TLegend(0.22,0.9,0.87,0.965)
leg.SetMargin(0.11)
leg.SetNColumns(1)

# ppb_stat040.Draw("Pz")
# ppb_syst040.Draw("P2")
# leg.AddEntry(ppb_syst040,"","pf")
# leg.AddEntry(pp_syst,"","pf")
leg.AddEntry(pbpb_syst,"","pf")


leg.SetEntrySeparation(0.2)


legT = ROOT.TLegend(0.64,0.17,0.98,0.42)
legT.SetMargin(0.14)
legT.SetBorderSize(0)
legT.AddEntry(s3_3body, s3_3body.GetTitle(), "LF")
legT.AddEntry(s3_2body, s3_2body.GetTitle(), "LF")
legT.AddEntry(s3_csm_1, s3_csm_1.GetTitle(), "LF")
legT.AddEntry(s3_csm_3, s3_csm_3.GetTitle(), "LF")
# legT.AddEntry(hp_ratio_csm_van)


legT2 = ROOT.TLegend(0.54,0.2,0.93, 0.31)
# legT2.SetMargin(0.14)
legT2.AddEntry(s3_csm_1, s3_csm_1.GetTitle(), "LF")
legT2.AddEntry(s3_csm_3, s3_csm_3.GetTitle(), "LF")
# legT.AddEntry(hp_ratio_csm_van)
leg.SetFillStyle(0)
legT.SetFillStyle(0)
legT2.SetFillStyle(0)



leg.Draw()
legT.Draw()
# legT2.Draw()





# leg.AddEntry(pbpb_syst,"","pf")
# leg.SetEntrySeparation(0.2)
# legT = ROOT.TLegend(0.18,0.87062 ,0.929085,0.985175)
# legT.SetNColumns(2)

# legT.SetEntrySeparation(0.2)
# legT.SetMargin(0.14)

# legT.AddEntry(s3_2body)
# legT.AddEntry(s3_3body)
# legT.AddEntry(s3_csm_1)
# legT.AddEntry(s3_csm_3)
# leg.SetFillStyle(0)
# legT.SetFillStyle(0)
# leg.Draw()
# legT.Draw()

cv.Draw()

cv.SaveAs("../Results/s3.pdf")
cv.SaveAs("../Results/s3.png")
cv.SaveAs("../Results/s3.eps")

