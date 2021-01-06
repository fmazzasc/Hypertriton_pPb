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
c0 = ROOT.RooRealVar('c0', 'constant c0', -2,2)
c1 = ROOT.RooRealVar('c1', 'constant c1', -2, 2)
hyp_mass = ROOT.RooRealVar('hyp_mass', 'hypertriton mass', 2.989, 2.993, 'GeV/c^{2}')
width = ROOT.RooRealVar('width', 'hypertriton width', 0., 0.008, 'GeV/c^{2}')
bkg_pdf_fit = ROOT.RooPolynomial('bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0, c1))
# define signal and background normalization
n1 = ROOT.RooRealVar('n1', 'n1 const', 0., 1, 'GeV')




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



    yield_gen_fit = ROOT.TH1D(f'yield_gen_fit_{eff}',"; (Generated - Fitted)/Generated signal", 40, -3, 3)

    bkg_pdf = bkg_file.Get(f'bkg_pdf/bkg_pdf_{eff}_040')
    cut_string = f"model_output>{cut}"
    len_data = len(np.array(df_data.query(cut_string + " and 2.96<m<3.04 ")["m"]))
    data = np.array(df.query(cut_string)["m"])
    roo_data = hp.ndarray2roo(data, mass)
    roo_hist = roo_data.binnedClone()

    signal_pdf = ROOT.RooHistPdf("histpdf1", "histpdf1", ROOT.RooArgSet(mass), roo_hist, 0)
    gen_function = ROOT.RooAddPdf('signal+bkg(toy)', 'signal+bkg(toy)', ROOT.RooArgList(signal_pdf, bkg_pdf), ROOT.RooArgList(pdf_frac))

    signal_pdf_fit = ROOT.RooHistPdf("histpdf_fit", "histpdf_fit", ROOT.RooArgSet(mass), roo_hist, 0)
    fit_function = ROOT.RooAddPdf('signal+bkg(fit)', 'signal+bkg(fit)', ROOT.RooArgList(signal_pdf_fit, bkg_pdf_fit), ROOT.RooArgList(n1))

    initial_signal = len_data*pdf_frac.getVal()
    print('INITIAL SIGNAL: ', initial_signal)


    workspace = ROOT.RooWorkspace()
    getattr(workspace, 'import')(gen_function)
    workspace.saveSnapshot(f'pdf_{eff}','c0')
    workspace.saveSnapshot(f'pdf_{eff}','c1')

    
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


        yield_gen_fit.Fill(float(initial_signal - fitted_signal)/initial_signal) 
        workspace.loadSnapshot(f'pdf_{eff}')
 
    yield_gen_fit.SetDirectory(0)

bkg_file.Close()






output_file = ROOT.TFile('../../Utils/toy_signal.root', 'recreate')
yield_gen_fit.Write()



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