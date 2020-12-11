import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mplhep as mpl

plt.style.use(mpl.style.ALICE)


abs_frac_list = ["", "_1.5", "_2", "_10"]
abs_array_list = []
for ind in abs_frac_list:
    abs_hist = uproot.open('../../Utils/AbsorptionHe3/absorption_pt/recPtHe3' + ind + ".root")['Reconstructed pT spectrum']
    abs_array_list.append(abs_hist.values)

pt_bins = abs_hist.edges


abs_perc_list = ['100%', '150%', '200%', '1000%']
fig = plt.figure()
for ind, perc in enumerate(abs_perc_list):
    plt.step(pt_bins, np.insert(abs_array_list[ind], 0, abs_array_list[ind][0]), label = perc)

plt.xlim((pt_bins[0],pt_bins[-1]))
legend = plt.legend(title=r'Percentage of $^{3}\overline{He}$ cross section', fontsize=14, loc='upper left')
legend.get_title().set_fontsize('16')
plt.xlabel(r'$\it{p}_{\mathrm{T}}$ $(GeV/\it{c})$')
plt.ylabel(r'$P_{abs}$')
plt.savefig('../../Results/abs_syst.png')


df = uproot.open('../../Tables/SignalTable_17d_mtexp.root')['SignalTable'].pandas.df()

mean_lrec = []
std_lrec = []


for pt_ind in range(len(pt_bins) - 1):
    pt_low = pt_bins[pt_ind]
    pt_up = pt_bins[pt_ind + 1]
    mean_lrec.append(np.mean(df.query('@pt_low<pt<@pt_up')['Lrec']))
    std_lrec.append(np.std(df.query('@pt_low<pt<@pt_up')['Lrec']))

plt.figure()
mean_lrec.insert(0, mean_lrec[0])
plt.step(pt_bins, mean_lrec)
plt.xlim((pt_bins[0],pt_bins[-1]))
plt.xlabel(r'$\it{p}_{\mathrm{T}}$ $(GeV/\it{c})$')
plt.ylabel('Average decay length')
plt.show()
