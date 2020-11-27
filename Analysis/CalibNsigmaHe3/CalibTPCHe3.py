import os
import sys
import pickle

import matplotlib
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import ROOT
import uproot
from ROOT import RooFit as rf
from scipy.interpolate import UnivariateSpline

sys.path.append("../")
import helpers as hp

matplotlib.use("pdf")
matplotlib.style.use(mplhep.style.ALICE)

ROOT.ROOT.EnableImplicitMT()
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

ROOT.gROOT.LoadMacro('RooGausDExp.cxx')
ROOT.gROOT.LoadMacro('CalibHe3.cc')

from ROOT import CalibHe3

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)



nsigma_calib_path = "../../Utils/CalibNsigmaHe3"
if not os.path.exists(nsigma_calib_path):
        os.makedirs(nsigma_calib_path)


mu_list = []
mu_error_list = []
sigma_list = []
sigma_error_list = []

df = uproot.open("../../Tables/DataTable_16qt_pass2.root")["DataTable"].pandas.df(["TPCsignalHe3", "TPCmomHe3","TPCnSigmaHe3", "V0CosPA", "NpidClustersHe3", "ArmenterosAlpha"])

mom_values = np.round(np.arange(0.55, 1.3, 0.05), 2)
output_file = ROOT.TFile(nsigma_calib_path + "/nsigmas.root",'recreate')

for he3_mom in mom_values:
    sel_df = df.query("NpidClustersHe3>155 and @he3_mom - 0.025<TPCmomHe3<@he3_mom + 0.025")
    bins = 60 if he3_mom < 1 else 30
    sup_edge = np.max(sel_df["TPCsignalHe3"])
    #### define low edge####
    if he3_mom < 1.1:
        res = np.histogram(sel_df["TPCsignalHe3"], bins=bins)
        counts = res[0]
        edges=res[1]
        bin_center = (edges[1:] + edges[:-1])/2
        bin_max = bin_center[counts==np.max(counts)]
        low_edge = bin_max + 50
    else:
        low_edge = 425

    tpc_histo = ROOT.TH1D(f"TPCsignalHe3_{he3_mom}", f"TPCsignalHe3_{he3_mom}", bins, low_edge, sup_edge)
    for signal in sel_df["TPCsignalHe3"]:
        tpc_histo.Fill(signal)

    # Define TPC variable
    ns = ROOT.RooRealVar(r'signalHe3',r'signalHe3',low_edge, sup_edge,'a. u.')

    # signal
    nSigTPC = ROOT.RooRealVar(r'N_{sig}','counts TPC signal',1e+2,0,1e+4)
    gaus_low = low_edge + 100 if he3_mom < 1 else low_edge
    sig_max = 200 if he3_mom < 1 else 40
    muSigTPC = ROOT.RooRealVar(r'#mu_{sig}','mu TPC signal', (low_edge + sup_edge)/2, gaus_low, sup_edge)
    sigmaSigTPC = ROOT.RooRealVar(r'#sigma_{sig}',' TPC signal', 10, sig_max)


    sigModelTPC = ROOT.RooGaussian('signal','signal function TPC', ns, muSigTPC, sigmaSigTPC)

    # bkg 1
    nBkg0TPC = ROOT.RooRealVar(r'N_{bkg0}','counts TPC bkg0', 0, 1e+4)
    tauBkg0TPC = ROOT.RooRealVar(r'#tau_{bkg0}','tau TPC bkg0',-10, 10)
    bkg0ModelTPC = ROOT.RooExponential('bkg0','bkg0 function TPC',ns,tauBkg0TPC)

    nBkg1TPC = ROOT.RooRealVar(r'N_{bkg1}','counts TPC bkg1', 0, 1e+4)
    tauBkg1TPC = ROOT.RooRealVar(r'#tau_{bkg1}','tau TPC bkg1',-10, 10)
    bkg1ModelTPC = ROOT.RooExponential('bkg1','bkg1 function TPC',ns,tauBkg1TPC)

    totModelTPC = ROOT.RooAddPdf('model','tot function TPC', ROOT.RooArgList(sigModelTPC, bkg0ModelTPC, bkg1ModelTPC),ROOT.RooArgList(nSigTPC,nBkg0TPC, nBkg1TPC))

    hTPChist = ROOT.RooDataHist(f"TPCsignalHe3_{he3_mom}", f"TPCsignalHe3_{he3_mom}", ROOT.RooArgList(ns), tpc_histo)

    frame, chi2,fitres = hp.plotData(ns, hTPChist, totModelTPC, 1.1, "SignalHe3", low_edge , sup_edge, f"SignalHe3_{he3_mom}")

    mu_list.append(muSigTPC.getVal())
    mu_error_list.append(muSigTPC.getError())
    sigma_list.append(sigmaSigTPC.getVal())
    sigma_error_list.append(sigmaSigTPC.getError())

    
    frame.Write()


mean_histo = ROOT.TH1D(f"Mean TPCsignal", f"Mean TPCsignal", len(mom_values), mom_values[0] - 0.025, mom_values[-1] + 0.025)
sigma_histo = ROOT.TH1D(f"Sigma TPCsignal", f"Sigma TPCsignal", len(mom_values), mom_values[0] - 0.025, mom_values[-1] + 0.025)

for iBin in range(1, mean_histo.GetNbinsX() + 1):
    mean_histo.SetBinContent(iBin, mu_list[iBin - 1])
    mean_histo.SetBinError(iBin, mu_error_list[iBin - 1])
    sigma_histo.SetBinContent(iBin, sigma_list[iBin - 1])
    sigma_histo.SetBinError(iBin, sigma_error_list[iBin - 1])

CalibHe3(mean_histo, sigma_histo, output_file)




