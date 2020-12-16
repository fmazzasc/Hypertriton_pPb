import sys
sys.path.append('../')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import mplhep
import ROOT

import helpers as hp

ROOT.gROOT.SetBatch()


pdf_fraction_list = np.load('../../Utils/pdf_fraction_040.npy')
toy_copies = 10000
# define working variable
mass = ROOT.RooRealVar('m', 'm_{^{3}He+#pi}', 2.96, 3.04, 'GeV/c^{2}')
# define signal parameters
c0 = ROOT.RooRealVar('c0', 'constant c0', 0., 100.)
hyp_mass = ROOT.RooRealVar('hyp_mass', 'hypertriton mass', 2.989, 2.993, 'GeV/c^{2}')
width = ROOT.RooRealVar('width', 'hypertriton width', 0., 0.008, 'GeV/c^{2}')
# define signal component
signal = ROOT.RooGaussian('signal', 'signal component pdf', mass, hyp_mass, width)
background = ROOT.RooPolynomial('bkg', 'pol0 bkg', mass, ROOT.RooArgList(c0))
# define signal and background normalization
n1 = ROOT.RooRealVar('n1', 'n1 const', 0., 1, 'GeV')
# define the fit funciton -> signal component + background component
fit_function = ROOT.RooAddPdf('model', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(n1))


df = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")
df_data = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_data.parquet.gzip")
score_cuts_array = np.load("../../Utils/Efficiencies/score_eff_syst_arr.npy")
bdt_eff_array = np.load("../../Utils/Efficiencies/bdt_eff_syst_arr.npy")
selected_bdt_eff = 0.72
bins = 500
mass.setBins(bins)

bkg_file = ROOT.TFile("../../Results/inv_mass_fits.root")

hist_list = []

toy_data_list = []
func_list = []


for pdf_frac, eff,cut in zip(pdf_fraction_list, bdt_eff_array, score_cuts_array):

    if eff!=0.72:
        continue


    n1.setVal(pdf_frac.getVal())


    # n1.setConstant(True)

    correction = ROOT.TH1D(f'Width_{eff}',"; Width (GeV/#it{c}^{2})", 80, 0.00, 0.008)
    yield_width = ROOT.TH2D(f'yield_width_{eff}',"; Width (GeV/#it{c}^{2}); (Generated - Fitted)/Generated signal", 80, 0.001, 0.008, 40, -3, 3)
    print(yield_width)
    bkg_pdf = bkg_file.Get(f'bkg_pdf/bkg_pdf_{eff}_040')


    cut_string = f"model_output>{cut}"
    len_data = len(np.array(df_data.query(cut_string + " and 2.96<m<3.04 ")["m"]))
    data = np.array(df.query(cut_string)["m"])
    roo_data = hp.ndarray2roo(data, mass)
    roo_hist = roo_data.binnedClone()
    signal_pdf = ROOT.RooHistPdf("histpdf1", "histpdf1", ROOT.RooArgSet(mass), roo_hist, 0)

    gen_function = ROOT.RooAddPdf('signal+bkg(toy)', 'signal+bkg(toy)', ROOT.RooArgList(signal_pdf, bkg_pdf), ROOT.RooArgList(pdf_frac))

    initial_signal = len_data*pdf_frac.getVal()
    print('INITIAL SIGNAL: ', initial_signal)


    workspace = ROOT.RooWorkspace()
    getattr(workspace, 'import')(gen_function)
    workspace.saveSnapshot(f'pdf_{eff}','c0')

    
    for toy in range(toy_copies):
        print(toy)
        wpdf = workspace.pdf(gen_function.GetName())
        toy_data = wpdf.generate(ROOT.RooArgSet(mass), len_data)
        
        fit_function.fitTo(toy_data)
        frame = mass.frame(36)
        toy_data.plotOn(frame)
        fit_function.plotOn(frame)
        fitted_signal = len_data*n1.getVal()
        

        if toy==15:

            toy_data_list.append(toy_data)
            func_list.append(fit_function)


        if width.getVal()<0.0077:
            correction.Fill(width.getVal())
            yield_width.Fill(width.getVal(), float(initial_signal - fitted_signal)/initial_signal)   
        workspace.loadSnapshot(f'pdf_{eff}')
 
    correction.SetDirectory(0)
    yield_width.SetDirectory(0)

bkg_file.Close()



Fcn1 = ROOT.TLine(0.00151, correction.GetMinimum(), 0.00151, correction.GetMaximum() + 4)
Fcn1.SetLineWidth(1)
Fcn1.SetLineStyle(2)
Fcn1.SetLineColor(ROOT.kRed - 3)


boxFcn1 = ROOT.TBox(0.00117, correction.GetMinimum(), 0.00185, correction.GetMaximum() + 4)
boxFcn1.SetFillColorAlpha(ROOT.kRed - 3, 0.35)
boxFcn1.SetFillStyle(3004)
boxFcn1.SetLineWidth(1)
boxFcn1.SetLineStyle(2)
boxFcn1.SetLineColor(ROOT.kRed - 3)


output_file = ROOT.TFile('../../Utils/width.root', 'recreate')
yield_width.Write()
cv = ROOT.TCanvas(f"width_0.72")
correction.SetMaximum(correction.GetMaximum() + 4)
correction.Fit('gaus', '', '', 0.002 - 0.001 , 0.002 + 0.001)
correction.Draw()
boxFcn1.Draw("same")
Fcn1.Draw('lsame')
leg1 = ROOT.TLegend(.68,0.26,.78,0.43)
leg1.SetFillStyle(0)
leg1.SetMargin(0.2)
leg1.SetBorderSize(0)
leg1.SetTextFont(42)
leg1.SetTextSize(0.025)
leg1.AddEntry(correction, "Width from 1000 toys", "l")
leg1.AddEntry(boxFcn1, "Width from data", "fl")
leg1.Draw('same')
cv.Write()


for toy,func in zip(toy_data_list, func_list):

    frame = mass.frame(36)
    frame.SetTitle("")

    toy.plotOn(frame)
    func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))
    # func.paramOn(frame)
    
    cv = ROOT.TCanvas(f"toy_example")
    frame.Draw()
    cv.Write()

output_file.Close()