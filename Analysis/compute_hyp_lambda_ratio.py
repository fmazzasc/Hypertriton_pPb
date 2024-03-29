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
import math
import ROOT
ROOT.gROOT.SetBatch()



##COMPUTE PRESELECTION-EFFICIENCY
df_rec = uproot.open("../Tables/SignalTable_20l2_mtexp.root")["SignalTable"].arrays(library="pd")
df_sim = uproot.open("../Tables/SignalTable_20l2_mtexp.root")["GenTable"].arrays(library="pd")

presel_eff = len(df_rec)/len(df_sim.query("-0.5<rapidity<0.5"))
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
    data040 = np.array(df.query(cut_string + " and 2.96<m<3.04 and centrality<=40 and abs(fZ) < 10")["m"])
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


###COMPUTE YIELD############################


corrected_counts = signal_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
corrected_error = error_array040/n_events040/branching_ratio/presel_eff/bdt_eff_array
Yield = np.float64(corrected_counts[bdt_eff_array==selected_bdt_eff]  / 2)
Yield = Yield/0.97   #absorption correction

stat_error = np.float64(corrected_error[bdt_eff_array==selected_bdt_eff] / 2 / 0.97)
# syst_error = float(np.std(corrected_counts) / corrected_counts[bdt_eff_array==selected_bdt_eff])
syst_error = 0.14*Yield
pt_shape_syst = 0.07*Yield
abs_syst = 0.041*Yield
br_syst = 0.09*Yield
syst_error = np.sqrt(syst_error**2 + pt_shape_syst**2 + abs_syst**2 + br_syst**2)
print("-------------------------------------")
sig = signal_array040[bdt_eff_array==selected_bdt_eff][0]
err_sig = error_array040[bdt_eff_array==selected_bdt_eff][0]
print(f"Signal = {sig:.3e} +- {err_sig:.3e} (stat.)")
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



lambdaAv040 = hp.computeAverage(lambdaVals,40)


# d$N$/d$\eta$ obtained by simple weighted average of the values published in https://arxiv.org/pdf/1910.14401.pdf
x_pPb040=np.array([29.4], dtype=np.float64)
xe_pPb040=np.array([0.6], dtype=np.float64)


hp_ratio_040 = np.array([2*Yield / lambdaAv040[0]], dtype=np.float64)
hp_ratiostat040 = np.array([hp_ratio_040[0] * hp.myHypot(stat_error / Yield, lambdaAv040[1] / lambdaAv040[0])], dtype=np.float64)
hp_ratiosyst040 = np.array([hp_ratio_040[0] * hp.myHypot(syst_error / Yield, lambdaAv040[2] / lambdaAv040[0])], dtype=np.float64)

print("-------------------------------------")
print(f"Hyp/Lambda 0-40% = {hp_ratio_040[0]:.2e} +- {hp_ratiostat040[0]:.2e} (stat.) +- {hp_ratiosyst040[0]:.2e} (syst.)")

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

hp_ratio_csm_3 = ROOT.TGraphErrors("../Utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.3dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")
hp_ratio_csm_1 = ROOT.TGraphErrors("../Utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")


hp_ratio_csm_1.SetLineColor(922)
hp_ratio_csm_1.SetLineWidth(2)
hp_ratio_csm_1.SetTitle("SHM, #it{Vc} = d#it{V}/d#it{y}")
hp_ratio_csm_1.SetMarkerSize(0)


hp_ratio_csm_3.SetLineColor(922)
hp_ratio_csm_3.SetLineWidth(2)
hp_ratio_csm_3.SetLineStyle(2)
hp_ratio_csm_3.SetMarkerSize(0)
hp_ratio_csm_3.SetTitle("SHM, #it{Vc} = 3d#it{V}/d#it{y}")

n = hp_ratio_csm_1.GetN()
grshade = ROOT.TGraph(2*n)
for i in range(n) : 
   grshade.SetPoint(i, hp_ratio_csm_3.GetPointX(i), hp_ratio_csm_3.GetPointY(i))
   grshade.SetPoint(n + i, hp_ratio_csm_1.GetPointX(n - i -1), hp_ratio_csm_1.GetPointY(n - i - 1))
   
grshade.SetFillColorAlpha(16, 0.571)
# grshade.SetFillStyle(3013)




hp_2body = ROOT.TGraphErrors("../Utils/ProdModels/coalescence/hp_ratio_2body_coal.csv","%lg %lg %lg")
hp_2body.SetLineColor(kBlueC)
hp_2body.SetMarkerColor(kBlueC)
hp_2body.SetTitle("2-body coalescence")
hp_2body.SetFillStyle(3145)
hp_2body.SetMarkerSize(0)
hp_2body.SetLineWidth(2)

hp_3body = ROOT.TGraphErrors("../Utils/ProdModels/coalescence/hp_ratio_3body_coal.csv","%lg %lg %lg")
hp_3body.SetLineColor(kAzureC)
hp_3body.SetMarkerColor(kAzureC)
hp_3body.SetTitle("3-body coalescence")
hp_3body.SetFillStyle(3014)
hp_3body.SetMarkerSize(0)
hp_3body.SetLineWidth(2)


hp_3body.SetFillColorAlpha(kAzureC, 0.571)
hp_2body.SetFillColorAlpha(kBlueC, 0.571)

mg = ROOT.TMultiGraph()
mg.Add(hp_2body)
mg.Add(hp_3body)



cv = ROOT.TCanvas("cv", "cv", 700,700)
cv.SetBottomMargin(0.145)
cv.SetLeftMargin(0.17)
cv.SetTopMargin(0.01)
cv.SetRightMargin(0.01)
frame=cv.DrawFrame(5., 1e-7, 200, 6e-6,";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}; {}_{#Lambda}^{3}H/p")

cv.SetLogx()
cv.SetLogy()
mg.Draw("4al same")
# grshade.Draw("f same")
hp_ratio_csm_1.Draw("L same")
hp_ratio_csm_3.Draw("L same")
# hp_ratio_csm_van.Draw("L same")
mg.GetYaxis().SetRangeUser(2e-7, 2e-5)
mg.GetXaxis().SetRangeUser(5, 3e3)
mg.GetXaxis().SetTitle('#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}')
mg.GetYaxis().SetTitle('{}_{#Lambda}^{3}H/#Lambda')
mg.GetXaxis().SetTitleOffset(1.1)
mg.GetYaxis().SetTitleOffset(1.16)
mg.GetYaxis().SetTitleSize(0.06)
mg.GetXaxis().SetTitleSize(0.06)
mg.GetYaxis().SetLabelSize(0.04)
mg.GetXaxis().SetLabelSize(0.04)



zero = np.array([0], dtype=np.float64)




ppb_stat040 = ROOT.TGraphErrors(1,x_pPb040,hp_ratio_040,zero,hp_ratiostat040)
ppb_stat040.SetLineColor(kRedC)
ppb_stat040.SetMarkerColor(kRedC)
ppb_stat040.SetMarkerStyle(20)
ppb_stat040.SetMarkerSize(1)
ppb_stat040.SetLineWidth(1)

ppb_syst040 = ROOT.TGraphErrors(1,x_pPb040, hp_ratio_040, xe_pPb040, hp_ratiosyst040)
ppb_syst040.SetTitle("ALICE p#font[122]{-}Pb, 0#font[122]{-}40%, #sqrt{#it{s}_{NN}} = 5.02 TeV")
ppb_syst040.SetLineColor(kRedC)
ppb_syst040.SetMarkerColor(kRedC)
ppb_syst040.SetFillStyle(0)
ppb_syst040.SetMarkerStyle(20)
ppb_syst040.SetLineWidth(1)

ppb_syst040.SetMarkerSize(1)



x = np.array([30.81], dtype=np.float64)
ex = np.array([0.44], dtype=np.float64)
y = np.array([5.49*10**(-7)], dtype=np.float64)
ey = np.array([1.34*10**(-7)], dtype=np.float64)
eys = np.array([0.72e-7], dtype=np.float64)
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
y = np.array([4*1.64e-6], dtype=np.float64)
ey = np.array([4*3.28e-7], dtype=np.float64)
eys = np.array([4*1.69e-7], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_stat.SetLineColor(ROOT.kBlue + 3)
pbpb_stat.SetMarkerColor(ROOT.kBlue +3)
pbpb_stat.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_stat.SetMarkerSize(1.5)
pbpb_stat.SetLineWidth(1)
pbpb_stat.Draw("Pz")


pbpb_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_syst.SetTitle("ALICE Pb#font[122]{-}Pb, 0#font[122]{-}10%, #sqrt{#it{s}_{NN}} = 2.76 TeV")
pbpb_syst.SetLineColor(ROOT.kBlue +3)
pbpb_syst.SetMarkerColor(ROOT.kBlue +3)
pbpb_syst.SetFillStyle(0)
pbpb_syst.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_syst.SetMarkerSize(1.5)
pbpb_syst.SetLineWidth(1)
pbpb_syst.Draw("P2")


leg = ROOT.TLegend(0.221,0.9,0.87,0.96)
leg.SetMargin(0.12)
leg.SetNColumns(1)

# ppb_stat040.Draw("Pz")
# ppb_syst040.Draw("P2")
# leg.AddEntry(ppb_syst040,"","pf")
# leg.AddEntry(pp_syst,"","pf")
leg.AddEntry(pbpb_syst,"","pf")




leg.SetEntrySeparation(0.2)
legT = ROOT.TLegend(0.62,0.17,0.97,0.44)
legT.SetMargin(0.14)
legT.SetBorderSize(0)
legT.AddEntry(hp_3body, hp_3body.GetTitle(), "LF")
legT.AddEntry(hp_2body,  hp_2body.GetTitle(), "LF")
legT.AddEntry(hp_ratio_csm_1,  hp_ratio_csm_1.GetTitle(), "LF")
legT.AddEntry(hp_ratio_csm_3, hp_ratio_csm_3.GetTitle(), "LF")
# legT.AddEntry(hp_ratio_csm_van)
leg.SetFillStyle(0)
legT.SetFillStyle(0)
leg.Draw()
legT.Draw()

pinfo = ROOT.TPaveText(0.2208,0.84, 0.47, 0.89, 'NDC')
pinfo.SetBorderSize(0)
pinfo.SetFillStyle(0)
# pinfo.SetTextAlign(30+3)
pinfo.SetTextFont(42)
pinfo.AddText('B.R. = 0.25 #pm 0.02')
pinfo.Draw()


cv.Draw()

cv.SaveAs("../Results/hl_ratio.pdf")
cv.SaveAs("../Results/hl_ratio.png")
cv.SaveAs("../Results/hl_ratio.eps")

# Branching ratio plot


# Branching ratio plot
cvBR = ROOT.TCanvas("brPlot","BR Plot",700,800)
cvBR.SetTopMargin(0.06)
cvBR.SetLeftMargin(0.15)
cvBR.SetRightMargin(0.04)
frame = cvBR.DrawFrame(0.25 - 0.022, 5.e-8, 0.25 + 0.022, 2e-6, ";B.R. (^{3}_{#Lambda}H #rightarrow ^{3}He + #pi^{-}) ; {}_{#Lambda}^{3}H/#Lambda #times B.R.")
frame.GetYaxis().SetTitleOffset(1.4)
cvBR.SetLogy()

npoints = 20000
xBR = np.linspace(0.25 - 0.022, 0.25 + 0.022,npoints)
C2 = 6.57e-7
C2e = 1.39e-7
C3 = 3.15e-7
C3e = 8.57e-8

CSM1 = 2.58e-06
CSM3 = 4.59e-06
CSM = 0.5 * (CSM3 + CSM1)
CSMe = 0.5 * (CSM3 - CSM1)

grData = ROOT.TGraphErrors(2)
grData.SetPoint(0, 0., hp_ratio_040[0] * 0.25)
grData.SetPointError(0, 0., math.hypot(hp_ratiostat040[0], hp_ratiosyst040[0]) * 0.25)
grData.SetPoint(1, 1., hp_ratio_040[0] * 0.25)
grData.SetPointError(1, 0., math.hypot(hp_ratiostat040[0], hp_ratiosyst040[0]) * 0.25)
grData.SetLineColor(kRedC)
grData.SetLineWidth(2)
grData.SetFillColorAlpha(kRedC,0.571)
grData.SetTitle("ALICE p#font[122]{-}Pb, 0#font[122]{-}40%, #sqrt{#it{s}_{NN}} = 5.02 TeV")

grCSM = ROOT.TGraphErrors(npoints)
grCSM1 = ROOT.TGraph(npoints)
grCSM3 = ROOT.TGraph(npoints)
grC2 = ROOT.TGraphErrors(npoints)
grC3 = ROOT.TGraphErrors(npoints)
for i,j in enumerate(xBR) :
  grC2.SetPoint(i, j, C2 * j)
  grC2.SetPointError(i, 0, C2e * j)
  grC3.SetPoint(i, j, C3 * j)
  grC3.SetPointError(i, 0, C3e * j)
  grCSM1.SetPoint(i, j, CSM1 * j)
  grCSM3.SetPoint(i, j, CSM3 * j)
  grCSM.SetPoint(i, j, CSM * j)
  grCSM.SetPointError(i, 0, CSMe * j)

grshade = ROOT.TGraph(2*npoints)
for i in range(npoints) : 
   grshade.SetPoint(i, grCSM3.GetPointX(i), grCSM3.GetPointY(i))
   grshade.SetPoint(npoints + i, grCSM1.GetPointX(npoints - i -1), grCSM1.GetPointY(npoints - i - 1))
 
grshade.SetFillColorAlpha(16, 0.571)

grC2.SetLineColor(kBlueC)
grC2.SetMarkerColor(kBlueC)
grC2.SetMarkerSize(0)
grC2.SetFillColorAlpha(kBlueC, 0.571)
grC2.SetFillStyle(3145)
grC2.SetTitle("2-body coalescence")
# grC2.SetLineStyle(8)
grC2.SetLineWidth(2)
grC3.SetLineColor(kAzureC)
grC3.SetMarkerColor(kAzureC)
grC3.SetFillColorAlpha(kAzureC, 0.571)
# grC3.SetLineStyle(8)
grC3.SetFillStyle(3014)
grC3.SetTitle("3-body coalescence")
grC3.SetLineWidth(2)
grC3.SetMarkerSize(0)
grCSM.SetLineColor(kGreyC)
grCSM.SetFillColorAlpha(kGreyC, 0.571)
grCSM.SetLineColor(kGreyC)
grCSM1.SetLineWidth(2)
grCSM1.SetTitle("SHM, #it{Vc} = d#it{V}/d#it{y}")
grCSM1.SetLineStyle(1)
grCSM1.SetLineColor(kGreyC)

grCSM3.SetTitle("SHM, #it{Vc} = 3d#it{V}/d#it{y}")
grCSM3.SetLineColor(922)
grCSM3.SetLineWidth(2)
grCSM3.SetLineStyle(2)
grCSM3.SetMarkerSize(0)


frame.GetXaxis().SetDecimals()
frame.GetXaxis().SetNdivisions(505)
frame.GetXaxis().SetTitleSize(0.06)
frame.GetYaxis().SetTitleSize(0.06)
frame.GetXaxis().SetTitleOffset(.85)
frame.GetYaxis().SetTitleOffset(1.16)



grC2.Draw("CL3")
grC3.Draw("CL3")
grData.Draw("L3")
grCSM1.Draw("L same")
grCSM3.Draw("L same")

frame.GetXaxis().SetTitleSize(0.06)
frame.GetYaxis().SetTitleSize(0.06)

legData = ROOT.TLegend(0.2, 0.94 , 0.83, 0.99)
legData.SetMargin(0.064)
legData.SetFillStyle(0)
legData.AddEntry(grData, grData.GetTitle(), "LF")
legData.Draw()

legBR = ROOT.TLegend(0.2,0.83 , 1, 0.93)
legBR.SetNColumns(2)
legBR.SetMargin(0.1)
legBR.SetFillStyle(0)
legBR.AddEntry(grC2)
legBR.AddEntry(grCSM1, grCSM1.GetTitle(), "LF")
legBR.AddEntry(grC3)
legBR.AddEntry(grCSM3, grCSM3.GetTitle(), "LF")
legBR.Draw()

cvBR.SaveAs("../Results/plotBR.pdf")
cvBR.SaveAs("../Results/plotBR.png")
cvBR.SaveAs("../Results/plotBR.eps")



file = ROOT.TFile('../Results/hl_ratio.root', 'recreate')
cv.Write()
cvBR.Write()
file.Close()