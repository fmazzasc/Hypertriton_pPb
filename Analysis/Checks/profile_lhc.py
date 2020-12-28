import ROOT

#Open the rootfile and get the workspace from the exercise_0
fInput = ROOT.TFile("../../Utils/Workspace.root")
ws = fInput.Get("Workspace")
ws.Print()
ws.var("c0").setConstant(1)
ws.var("c1").setConstant(1)

model = ws.obj("ModelConfig")

sig_frac = ws.var("n1")
poi = ROOT.RooArgSet(sig_frac)
nullParams = poi.snapshot()
nullParams.setRealValue("n1",0.)

plc = ROOT.RooStats.ProfileLikelihoodCalculator()
plc.SetData(ws.data("data"))
plc.SetModel(model)
plc.SetParameters(poi)
plc.SetNullParameters(nullParams)

htr = plc.GetHypoTest()


#We get a HypoTestResult out of the calculator, and we can query it.
htr = plc.GetHypoTest()

print("-------------------------------------------------")
print("The p-value for the null is ", htr.NullPValue())
print("Corresponding to a signifcance of ", htr.Significance())
print("-------------------------------------------------")

#PyROOT sometimes fails cleaning memory, this helps
del plc