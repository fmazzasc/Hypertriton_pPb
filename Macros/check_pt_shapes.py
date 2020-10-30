import uproot
import numpy as np
import mplhep as mpl
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
import pandas as pd

plt.style.use(mpl.style.ALICE)
ptexp = uproot.open("../Tables/SignalTable_17d_ptexp.root")["SignalTable"].pandas.df(["pt"])
mtexp = uproot.open("../Tables/SignalTable_17d_mtexp.root")["SignalTable"].pandas.df(["pt"])
bol = uproot.open("../Tables/SignalTable_17d_bol.root")["SignalTable"].pandas.df(["pt"])
old_bw = pd.read_parquet("../Tables/SignalTable_old_bw.parquet", columns=["HypCandPt"])

fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(ptexp, bins=500, alpha=0.3, label=r"$\it{p}_{\mathrm{T}}$ exponential", histtype="step")
ax.hist(mtexp, bins=500, alpha=0.3, label=r"$\it{m}_{\mathrm{T}}$ exponential", histtype="step")
ax.hist(bol, bins=500, alpha=0.3, label= "Boltzmann", histtype="step")
ax.hist(old_bw, bins=500, alpha=0.3, label= "Pb-Pb cent reweighted BW", histtype="step")

handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
plt.legend(handles=new_handles, labels=labels)

plt.xlabel(r"$\it{p}_{\mathrm{T}}$ $(\mathrm{GeV}/\it{c})$")
plt.ylabel("Counts")
plt.xlim(0,10)
plt.savefig("../Results/pt_shapes.png")
plt.show()
