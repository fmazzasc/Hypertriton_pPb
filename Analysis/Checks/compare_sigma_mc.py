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


#df = uproot.open("../../Tables/SignalTable_17d_mtexp.root")['SignalTable'].pandas.df(['m'])
df = pd.read_parquet("../../Utils/ReducedDataFrames/selected_df_mc.parquet.gzip")

hist = ROOT.TH1D('after ML', 'after ML; Mass(GeV/c^{2}); Counts', 200, 2.975, 3.01)

for mass in df['m']:
    hist.Fill(mass)

ff = ROOT.TFile('width_after_ML.root', 'recreate')

hist.Fit('gaus', '', '', 2.98742, 2.992802)
hist.Write()