
import uproot
import numpy as np
import aghast
import ROOT

df_rec = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df()
df_sim = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["GenTable"].pandas.df()

df_rec2 = uproot.open("../../Tables/SignalTable_20l2_mtexp.root")["SignalTable"].pandas.df()
df_sim2 = uproot.open("../../Tables/SignalTable_20l2_mtexp.root")["GenTable"].pandas.df()


pt_bins = np.linspace(0.5, 10, 50)

rec_hist = np.histogram(df_rec['pt'], pt_bins)
sim_hist = np.histogram(df_sim['pt'], pt_bins)


rec_hist2 = np.histogram(df_rec2['pt'], pt_bins)
sim_hist2 = np.histogram(df_sim2['pt'], pt_bins)

ghost_hist = aghast.from_numpy(rec_hist)
rec_hist = aghast.to_root(ghost_hist, 'Efficiency_vs_pt')
rec_hist.GetXaxis().SetTitle("#it{p}_{T} GeV/c");
ghost_hist = aghast.from_numpy(sim_hist)
sim_hist = aghast.to_root(ghost_hist, 'sim_hist')
rec_hist.Divide(sim_hist)

ghost_hist2 = aghast.from_numpy(rec_hist2)
rec_hist2 = aghast.to_root(ghost_hist2, 'Efficiency_vs_pt_2')
rec_hist2.GetXaxis().SetTitle("#it{p}_{T} GeV/c");
ghost_hist2 = aghast.from_numpy(sim_hist2)
sim_hist2 = aghast.to_root(ghost_hist2, 'sim_hist2')
rec_hist2.Divide(sim_hist2)
rec_hist2.SetLineColor(ROOT.kRed)

leg = ROOT.TLegend(0.15,0.75,0.7,0.85)
leg.AddEntry(rec_hist)
leg.AddEntry(rec_hist2)

file = ROOT.TFile('eff_vs_pt.root', 'recreate')

cv = ROOT.TCanvas(f"eff")
rec_hist2.Draw()
rec_hist.Draw('same')
leg.Draw()

cv.Write()

file.Close()