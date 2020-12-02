import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mplhep as mpl
matplotlib.use("pdf")
plt.style.use(mpl.style.ALICE)

gen_ct = uproot.open("../Tables/SignalTable_17d_mtexp.root")["GenTable"].array(['ct'])
gen_ct = gen_ct[gen_ct < 35]
rand_arr = np.random.rand(len(gen_ct))

abs_frac_list = ["", "_1.5", "_2", "_10"]
abs_array_list = []
for ind in abs_frac_list:
    abs_hist = uproot.open('../Utils/AbsorptionHe3/absorption_ct/recCtHe3' + ind + ".root")['Reconstructed ct spectrum']
    abs_array_list.append(abs_hist.values)

ct_bins = abs_hist.edges


abs_perc_list = ['100%', '150%', '200%', '1000%']
fig = plt.figure()
for ind, perc in enumerate(abs_perc_list):
    plt.step(ct_bins, np.append(abs_array_list[ind], abs_array_list[ind][-1]), label = perc)

plt.xlim((ct_bins[0],ct_bins[-1]))
legend = plt.legend(title=r'Percentage of $^{3}\overline{He}$ cross section', fontsize=14, loc='upper right', bbox_to_anchor=(0.85, 0.94))
legend.get_title().set_fontsize('16')
plt.xlabel(r'$\it{c}t$ (cm)')
plt.ylabel(r'$P_{abs}$')
plt.savefig('../Results/abs_syst.png')


n_abs_array = np.zeros(4)

for i in range(len(ct_bins) - 1):
    ct_mask = np.logical_and(gen_ct > ct_bins[i], gen_ct <= ct_bins[i+1])
    ct_abs_prob = rand_arr[ct_mask]
    for j in range(len(n_abs_array)):
        abso = np.sum(ct_abs_prob < abs_array_list[j][i])
        gen = len(gen_ct[ct_mask])
        n_abs_array[j] += np.sum(ct_abs_prob < abs_array_list[j][i])

n_abs_array = n_abs_array / len(gen_ct)
print(n_abs_array)
print("Correction: ", round(n_abs_array[0], 4))
print("Absorption syst: ", round(100*(np.max(n_abs_array) - np.min(n_abs_array)), 2), "%")
