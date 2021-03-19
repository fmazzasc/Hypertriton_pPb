import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import mplhep
import ROOT
ROOT.gROOT.SetBatch()

def ndarray2roo(ndarray, var):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(
        'data', 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo

bins=1000


df_rec = uproot.open("~/Hypertriton_pPb/Tables/SignalTable_20l2_mtexp.root")["SignalTable"].arrays(library="pd")
mc_data = np.array(df_rec.query("2.96<m<3.04")["m"])[0:50000]
mass = ROOT.RooRealVar('m', 'm_{^{3}He+#pi}', 2.96, 3.04, 'GeV/c^{2}')
deltaMass = ROOT.RooRealVar("deltaM", '#Deltam', -0.06, 0.06, 'GeV/c^{2}')
shiftedMass = ROOT.RooAddition("mPrime", "m + #Deltam", ROOT.RooArgList(mass, deltaMass))
mc_data = ndarray2roo(mc_data, mass)


signal = ROOT.RooKeysPdf('signal', 'signal', shiftedMass, mass, mc_data, ROOT.RooKeysPdf.MirrorBoth, 2)

ffile = ROOT.TFile("c../../Results/heck_kde.root", "recreate")

frame = mass.frame(bins)
frame.SetTitle("")

mc_data.plotOn(frame)
signal.plotOn(frame)

frame.Write()

ffile.Close()