
import uproot
import numpy as np
import aghast
import ROOT

df_rec = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df()
df_sim = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["GenTable"].pandas.df()

print(len(df_rec)/len(df_sim))
print(len(df_sim)/len(df_sim))

pt_bins = np.linspace(0.5, 6, 50)

rec_hist = np.histogram(df_rec['pt'], pt_bins)
sim_hist = np.histogram(df_sim['pt'], pt_bins)


ghost_hist = aghast.from_numpy(rec_hist)
rec_hist = aghast.to_root(ghost_hist, 'Efficiency_vs_pt')
rec_hist.GetXaxis().SetTitle("#it{p}_{T} GeV/c");
ghost_hist = aghast.from_numpy(sim_hist)
sim_hist = aghast.to_root(ghost_hist, 'sim_hist')
rec_hist.Divide(sim_hist)

cv = ROOT.TCanvas(f"eff")
rec_hist.Draw()

cv.SaveAs('eff.png')