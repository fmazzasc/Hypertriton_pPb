import uproot
import numpy as np
import mplhep as mpl
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
matplotlib.use("pdf")

plt.style.use(mpl.style.ALICE)
ptexp = uproot.open("../../Tables/SignalTable_17d_ptexp.root")["SignalTable"].pandas.df(["pt"])
mtexp = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df(["pt"])
bol = uproot.open("../../Tables/SignalTable_17d_bol.root")["SignalTable"].pandas.df(["pt"])
bw = uproot.open("../../Tables/SignalTable_17d_bw.root")["SignalTable"].pandas.df(["pt"])

fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(ptexp, bins=500, label=r"$\it{p}_{\mathrm{T}}$ exponential", histtype="step", color="orange")
ax.hist(mtexp, bins=500, label=r"$\it{m}_{\mathrm{T}}$ exponential", histtype="step", color="red")
ax.hist(bol, bins=500, label= "Boltzmann", histtype="step", color="blue")
ax.hist(bw, bins=500, label= "Blast Wave", histtype="step")

handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
plt.legend(handles=new_handles, labels=labels)

plt.xlabel(r"$\it{p}_{\mathrm{T}}$ $(\mathrm{GeV}/\it{c})$")
plt.ylabel("Counts")
plt.xlim(0,10)
plt.savefig("../../Results/pt_shapes.png")

###Compute Syst due to pT shape
ptexp_gen = uproot.open("../../Tables/SignalTable_17d_ptexp.root")["GenTable"].pandas.df(["pt","rapidity"]).query("abs(rapidity)<0.5")["pt"]
mtexp_gen = uproot.open("../../Tables/SignalTable_17d_mtexp.root")["GenTable"].pandas.df(["pt","rapidity"]).query("abs(rapidity)<0.5")["pt"]
bol_gen = uproot.open("../../Tables/SignalTable_17d_bol.root")["GenTable"].pandas.df(["pt","rapidity"]).query("abs(rapidity)<0.5")["pt"]
bw_gen = uproot.open("../../Tables/SignalTable_17d_bw.root")["GenTable"].pandas.df(["pt","rapidity"]).query("abs(rapidity)<0.5")["pt"]


print(len(bw))
print(len(bw_gen))
eff = np.array([len(ptexp)/len(ptexp_gen), len(mtexp)/len(mtexp_gen), len(bol)/len(bol_gen), len(bw)/len(bw_gen)])
print(eff)

print("Maximum difference: ", np.max(eff) - np.min(eff))
print("pT shape syst: ", 100*(np.std(eff)/np.max(eff)), "%")

