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


df_orig = uproot.open("../../Tables/SignalTable_20g7_mtexp.root")['SignalTable'].pandas.df()
df_ml = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")

standard_selection = 'V0CosPA > 0.9999 and NpidClustersHe3 > 80 and He3ProngPt > 1.8 and pt > 2 and pt < 10 and PiProngPt > 0.15 and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and TPCnSigmaHe3 < 3.5 and TPCnSigmaHe3 > -3.5 and ProngsDCA < 1'
df_std = df_orig.query(standard_selection)

hist_orig = ROOT.TH1D('Raw', '; Mass(GeV/c^{2}); Counts', 200, 2.975, 3.01)
hist_ml = ROOT.TH1D('After ML', '; Mass(GeV/c^{2}); Counts', 200, 2.975, 3.01)
hist_std = ROOT.TH1D('After linear selections', '; Mass(GeV/c^{2}); Counts', 200, 2.975, 3.01)

for mass in df_orig['m']:
    hist_orig.Fill(mass)

for mass in df_ml['m']:
    hist_ml.Fill(mass)

for mass in df_std['m']:
    hist_std.Fill(mass)

hist_orig.Fit('gaus', '', '', 2.99109 - 0.002, 2.99109 + 0.002)
hist_ml.Fit('gaus', '', '', 2.98742, 2.992802)
hist_std.Fit('gaus', '', '', 2.99109 - 0.002, 2.99109 + 0.002)



ff = ROOT.TFile('../../Results/width_comparison.root', 'recreate')
hist_orig.Write()
hist_std.Write()
hist_ml.Write()

