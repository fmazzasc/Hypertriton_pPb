import os
import ROOT
from ROOT import RooFit as rf
import numpy as np
import matplotlib.pyplot as plt
import pickle

kBlueC  = ROOT.TColor.GetColor("#2077b4");
kRedC  = ROOT.TColor.GetColor("#d62827");
kGreenC  = ROOT.TColor.GetColor("#2ba02b");
kOrangeC  = ROOT.TColor.GetColor("#ff7f00");
kVioletC  = ROOT.TColor.GetColor("#9467bd");
kPinkC  = ROOT.TColor.GetColor("#e377c1");
kGreyC  = ROOT.TColor.GetColor("#7f7f7f");
kBrownC  = ROOT.TColor.GetColor("#8c564c");
kAzureC  = ROOT.TColor.GetColor("#18becf");
kGreenBC  = ROOT.TColor.GetColor("#1f78b4");
kBlueC  = ROOT.TColor.GetColor("#1f78b4");


def array2hist(counts, hist):
    for iBin in range(1, hist.GetNbinsX() + 1):
        hist.SetBinContent(iBin, counts[iBin-1])
        hist.SetBinError(iBin, np.sqrt(counts[iBin-1]))


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


def significance_error(signal, background, signal_error=None, background_error=None):

    if not signal_error:
        signal_error = np.sqrt(signal + 1e-10)
    if not background_error:
        background_error = np.sqrt(background + 1e-10)
    sb = signal + background + 1e-10
    sb_sqrt = np.sqrt(sb)

    s_propag = (sb_sqrt + signal / (2 * sb_sqrt))/sb * signal_error
    b_propag = signal / (2 * sb_sqrt)/sb * background_error

    if signal+background == 0:
        return 0

    return np.sqrt(s_propag * s_propag + b_propag * b_propag)


def unbinned_mass_fit(data, eff, bkg_model, output_dir, bkg_dir, cent_class, pt_range, ct_range, split, cent_string = '', bins=38, ws_name=''):

    # define working variable
    mass = ROOT.RooRealVar('m', 'm_{^{3}He+#pi}', 2.96, 3.04, 'GeV/c^{2}')

    # define signal parameters
    hyp_mass = ROOT.RooRealVar(
        'hyp_mass', 'hypertriton mass', 2.96, 3.04, 'GeV/c^{2}')
    width = ROOT.RooRealVar('width', 'hypertriton width',
                            0.001, 0.003, 'GeV/c^{2}')

    # define signal component
    signal = ROOT.RooGaussian(
        'signal', 'signal component pdf', mass, hyp_mass, width)

    # define background parameters
    slope = ROOT.RooRealVar('slope', 'exponential slope', -100., 100.)

    c0 = ROOT.RooRealVar('c0', 'constant c0', -2,2)
    c1 = ROOT.RooRealVar('c1', 'constant c1', -2, 2)
    c2 = ROOT.RooRealVar('c2', 'constant c2', -2, 2)

    # define background component depending on background model required
    if bkg_model == 'pol0':
        background = ROOT.RooPolynomial(
            'bkg', 'pol0 bkg', mass, ROOT.RooArgList(c0))

    if bkg_model == 'pol1':
        background = ROOT.RooPolynomial(
            'bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0, c1))

    if bkg_model == 'pol2':
        background = ROOT.RooPolynomial(
            'bkg', 'pol2 for bkg', mass, ROOT.RooArgList(c0, c1, c2))

    if bkg_model == 'expo':
        background = ROOT.RooExponential('bkg', 'expo for bkg', mass, slope)

    # define signal and background normalization
    n_sig = ROOT.RooRealVar('n_sig', 'n_sig', 0., 1000, 'GeV')
    n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', 0., 1000, 'GeV')



    # define the fit funciton -> signal component + background component
    fit_function = ROOT.RooAddPdf('model', 'signal + background',
                                  ROOT.RooArgList(signal, background), ROOT.RooArgList(n_sig, n_bkg))

    # convert data to RooData
    roo_data = ndarray2roo(data, mass)

    # fit data
    fit_function.fitTo(roo_data, rf.Range(
        2.96, 3.04))

    # plot the fit
    frame = mass.frame(bins)
    frame.SetTitle("")
    frame.GetXaxis().SetNdivisions(510)

    roo_data.plotOn(frame)
    fit_function.plotOn(frame, rf.LineColor(kBlueC))
    #fit_function.plotOn(frame, rf.Components('signal'), rf.LineStyle(ROOT.kDotted), rf.LineColor(ROOT.kRed))
    fit_function.plotOn(frame, rf.Components(
        'bkg'), rf.LineStyle(9), rf.LineColor(kOrangeC))

    # add info to plot
    nsigma = 3
    mu = hyp_mass.getVal()
    mu_error = hyp_mass.getError()
    sigma = width.getVal()
    sigma_error = width.getError()

    # compute significance
    mass.setRange('signal region',  mu -
                  (nsigma * sigma), mu + (nsigma * sigma))
    
    range_lower = mu - (nsigma * sigma)
    range_upper = mu + (nsigma * sigma)
    
    n_sig_val = n_sig.getVal()
    n_bkg_val = n_bkg.getVal()
    


    signal_counts = signal.createIntegral(ROOT.RooArgSet(
        mass), ROOT.RooArgSet(mass), 'signal region').getVal() * n_sig_val
    signal_error = signal.createIntegral(ROOT.RooArgSet(
        mass), ROOT.RooArgSet(mass), 'signal region').getVal() * (n_sig.getError())
    

    background_counts = background.createIntegral(ROOT.RooArgSet(
        mass), ROOT.RooArgSet(mass), 'signal region').getVal() * n_bkg_val
    background_error = background.createIntegral(ROOT.RooArgSet(
        mass), ROOT.RooArgSet(mass), 'signal region').getVal() * (n_bkg.getError())


    signif = signal_counts / np.sqrt(signal_counts + background_counts + 1e-10)
    signif_error = significance_error(signal_counts, background_counts, signal_error, background_error)

    pinfo = ROOT.TPaveText(0.537, 0.674, 0.837, 0.775, 'NDC')
    pinfo.SetBorderSize(0)
    pinfo.SetFillStyle(0)
    pinfo.SetTextAlign(30+3)
    pinfo.SetTextFont(42)
    # pinfo.SetTextSize(12)

    decay_label = {
        '': '{}^{3}_{#Lambda}H#rightarrow ^{3}He #pi^{-} + c.c.',
        '_matter': '{}^{3}_{#Lambda}H#rightarrow ^{3}He #pi^{-}',
        '_antimatter': '{}^{3}_{#bar{#Lambda}}#bar{H}#rightarrow ^{3}#bar{He}#pi^{+}',
    }

    string_list = []

    string_list.append(f'ALICE Internal, p-Pb 2016 + 2013, {cent_class[0]}-{cent_class[1]}%')
    string_list.append(
        decay_label[split] + '     %i #leq #it{p}_{T} < %i GeV/#it{c} ' % (pt_range[0], pt_range[1]))
    string_list.append(f'#mu = {mu*1000:.2f} #pm {mu_error*1000:.2f} ' + 'MeV/#it{c}^{2}')
    string_list.append(f'#sigma = {sigma*1000:.2f} #pm {sigma_error*1000:.2f} ' + 'MeV/#it{c}^{2}')

    if roo_data.sumEntries() > 0:
        string_list.append('#chi^{2} / NDF = ' + f'{frame.chiSquare(4):.2f}')

    string_list.append(f'Significance ({nsigma:.0f}#sigma) = {signif:.1f} #pm {signif_error:.1f}')
    string_list.append(f'S ({nsigma:.0f}#sigma) = {signal_counts:.1f} #pm {signal_error:.1f}')
    string_list.append(f'B ({nsigma:.0f}#sigma) = {background_counts:.1f} #pm {background_error:.1f}')

    if background_counts > 0:

        ratio = signal_counts / background_counts
        string_list.append(f'S/B ({nsigma:.0f}#sigma) = {ratio:.2f}')

    for s in string_list:
        pinfo.AddText(s)

    frame.addObject(pinfo)
    if output_dir != '':
        output_dir.cd()
    binning = 1000*((3.04-2.96)/bins)
    stry= f"Events/({binning:.2f}"
    stry += "MeV/#it{c}^{2})"
    frame.SetYTitle(stry)
    cv = ROOT.TCanvas(f"cv_gaus_{round(eff,2)}_{cent_string}")
    frame.Draw()
    if output_dir != '':
        cv.Write()
    
    if ws_name != '':
        w = ROOT.RooWorkspace(ws_name)
        mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
        mc.SetPdf(fit_function)
        mc.SetParametersOfInterest(ROOT.RooArgSet(n_sig))
        mc.SetObservables(ROOT.RooArgSet(mass))
        if bkg_model=='pol1':
            w.defineSet("nuisParams","n_bkg,c0,c1")
            mc.SetNuisanceParameters(w.set('nuisParams'))
        if bkg_model=='expo':
            w.defineSet("nuisParams","slope")
            mc.SetNuisanceParameters(w.set('nuisParams'))
        getattr(w,'import')(mc)
        getattr(w,'import')(roo_data)
        if not os.path.exists('../Utils/Workspaces'):
            os.makedirs('../Utils/Workspaces')
        w.writeToFile(f'../Utils/Workspaces/{ws_name}.root', True)


    # bkg_dir.cd()
    # background.SetName(f"bkg_pdf_{round(eff,2)}_{cent_string}")
    # background.SetTitle(f"bkg_pdf_{round(eff,2)}_{cent_string}")
    # background.Write()

    return signal_counts, signal_error, signif, signif_error, mu, mu_error, sigma, sigma_error, range_lower, range_upper


def unbinned_mass_fit_mc(data, eff, bkg_model, mc_data, output_dir, bkg_dir, cent_class, pt_range, ct_range, split, cent_string = '', bins=38, ws_name=''):

    # define working variable
    mass = ROOT.RooRealVar('m', 'm_{^{3}He+#pi}', 2.96, 3.04, 'GeV/c^{2}')
    # mass.setVal(3.1)

    deltaMass = ROOT.RooRealVar("deltaM", '#Deltam', -0.06, 0.06, 'GeV/c^{2}')
    shiftedMass = ROOT.RooAddition("mPrime", "m + #Deltam", ROOT.RooArgList(mass, deltaMass))

    # define signal component
    #hist_signal = ROOT.RooDataHist('signal_hist', 'signal_hist', ROOT.RooArgList(mass), signal_hist)
    # signal = ROOT.RooHistPdf("histpdf1", "histpdf1", ROOT.RooArgList(shiftedMass), ROOT.RooArgList(mass), hist_signal, 0)
    mc_data = ndarray2roo(mc_data, mass)

    signal = ROOT.RooKeysPdf('signal', 'signal', shiftedMass, mass, mc_data, ROOT.RooKeysPdf.MirrorBoth, 2)

    slope = ROOT.RooRealVar('slope', 'exponential slope', -100., 100.)
    c0 = ROOT.RooRealVar('c0', 'constant c0', -2,2)
    c1 = ROOT.RooRealVar('c1', 'constant c1', -2, 2)
    c2 = ROOT.RooRealVar('c2', 'constant c2', -2, 2)

    # define background component depending on background model required
    if bkg_model == 'pol0':
        background = ROOT.RooPolynomial(
            'bkg', 'pol0 bkg', mass, ROOT.RooArgList())

    if bkg_model == 'pol1':
        background = ROOT.RooPolynomial(
            'bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0, c1))

    if bkg_model == 'pol2':
        background = ROOT.RooPolynomial(
            'bkg', 'pol2 for bkg', mass, ROOT.RooArgList(c0, c1, c2))

    if bkg_model == 'expo':
        background = ROOT.RooExponential('bkg', 'expo for bkg', mass, slope)


    n_sig = ROOT.RooRealVar('n_sig', 'n_sig', 0., 1000, 'GeV')
    n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', 0., 1000, 'GeV')



    # define the fit funciton -> signal component + background component
    fit_function = ROOT.RooAddPdf('model', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(n_sig, n_bkg))    

    # convert data to RooData
    roo_data = ndarray2roo(data, mass)

    # fit data
    fit_function.fitTo(roo_data, rf.Range(
        2.96, 3.04))

    # plot the fit
    frame = mass.frame(bins)
    frame.SetTitle("")
    frame.GetYaxis().SetTitleOffset(0.5)
    frame.GetYaxis().SetTitleSize(0.07)


    roo_data.plotOn(frame, rf.MarkerSize(1.5))
    fit_function.plotOn(frame, rf.Components(
        'bkg'), rf.LineStyle(9), rf.LineColor(kOrangeC), rf.LineWidth(3))
    fit_function.plotOn(frame, rf.LineColor(kBlueC), rf.LineWidth(3), rf.DrawOption("same"))
    
    #fit_function.plotOn(frame, rf.Components('signal'), rf.LineStyle(ROOT.kDotted), rf.LineColor(ROOT.kRed))


    # add info to plot
    
    signal_counts = n_sig.getVal()
    signal_error = n_sig.getError()
    background_counts = n_bkg.getVal()
    background_error = n_bkg.getError()


    pinfo = ROOT.TPaveText(0.537, 0.574, 0.8, 0.875, 'NDC')
    pinfo.SetBorderSize(0)
    pinfo.SetFillStyle(0)
    pinfo.SetTextAlign(30+3)
    pinfo.SetTextFont(42)


    decay_label = {
        '': '{}^{3}_{#Lambda}H#rightarrow ^{3}He #pi^{-} + c.c.',
        '_matter': '{}^{3}_{#Lambda}H#rightarrow ^{3}He #pi^{-}',
        '_antimatter': '{}^{3}_{#bar{#Lambda}}#bar{H}#rightarrow ^{3}#bar{He}#pi^{+}',
    }

    string_list = []
    string_list.append("ALICE Internal p-Pb, 0-40%, #sqrt{#it{s}_{NN}}=5.02 TeV")   
    string_list.append(f'Signal = {signal_counts:.1f} #pm {signal_error:.1f}')   
    # string_list.append(f'Significance ({3:.0f}#sigma) = {signif:.1f} #pm {signif_error:.1f}')


    if background_counts > 0:

        ratio = signal_counts / background_counts
        # string_list.append(f'S/B ({3:.0f}#sigma) = {ratio:.2f}')

    for s in string_list:
        pinfo.AddText(s)

    # frame.addObject(pinfo)
    if output_dir != '':
        output_dir.cd()
    binning = 1000*((3.04-2.96)/bins)
    stry= f"Events/({binning:.2f}"
    stry += "MeV/#it{c}^{2})"
    frame.SetYTitle(stry)
    cv = ROOT.TCanvas(f"cv_templ_{round(eff,2)}_{bkg_model}_{cent_string}")
    frame.SetName(f"frame_{eff}")
    frame.Write()
    frame.Draw()

    if output_dir != '':
        cv.Write()
    if bkg_dir != '':
        bkg_dir.cd()
        background.SetName(f"bkg_pdf_{round(eff,2)}_{cent_string}")
        background.SetTitle(f"bkg_pdf_{round(eff,2)}_{cent_string}")
        background.Write()

    # if ws_name != '':
    #     w = ROOT.RooWorkspace(ws_name)
    #     mc = ROOT.RooStats.ModelConfig("ModelConfig",w)
    #     mc.SetPdf(fit_function)
    #     mc.SetParametersOfInterest(ROOT.RooArgSet(n_sig))
    #     mc.SetObservables(ROOT.RooArgSet(mass))
    #     if bkg_model=='pol1':
    #         w.defineSet("nuisParams","n_bkg,c0,c1")
    #         mc.SetNuisanceParameters(w.set('nuisParams'))
    #     if bkg_model=='expo':
    #         w.defineSet("nuisParams","slope")
    #         mc.SetNuisanceParameters(w.set('nuisParams'))
    #     getattr(w,'import')(mc)
    #     getattr(w,'import')(roo_data)
    #     if not os.path.exists('../Utils/Workspaces'):
    #         os.makedirs('../Utils/Workspaces')
    #     w.writeToFile(f'../Utils/Workspaces/{ws_name}.root', True)
    
    print("mass:", mass.getVal())
    print("shifted_mass:", shiftedMass.getVal())
    print("delta_mass:", deltaMass.getVal())


    return signal_counts, signal_error, deltaMass.getVal(), deltaMass.getError()


def fit_hist(
        histo, cent_class, pt_range, ct_range, nsigma=3, model="pol0", fixsigma=-1, sigma_limits=None, mode=2, split=''):
    # canvas for plotting the invariant mass distribution
    cv = ROOT.TCanvas(f'cv_{histo.GetName()}')

    # define the number of parameters depending on the bkg model
    if isinstance(model, list):
        n_bkgpars = len(model)
        fit_tpl = ROOT.TF1('fitTpl', f'pol2(0) + gausn(3) + gausn({n_bkgpars})', 0, 5)
        bkg_tpl = ROOT.TF1('bkgTpl', f'pol2(0) + gausn(3)', 0, 5)
        for i in range(n_bkgpars):
            fit_tpl.SetParName(i, f'B_{i}')
            fit_tpl.FixParameter(i, model[i])
        

    
    else:
        if 'pol' in str(model):
            n_bkgpars = int(model[3]) + 1
        elif 'expo' in str(model):
            n_bkgpars = 2
        else:
            print(f'Unsupported model {model}')

        fit_tpl = ROOT.TF1('fitTpl', f'{model}(0)+gausn({n_bkgpars})', 0, 5)
        bkg_tpl = ROOT.TF1('fitTpl', f'{model}(0)', 0, 5)

    # redefine parameter names for the bkg_model
        for i in range(n_bkgpars):
            fit_tpl.SetParName(i, f'B_{i}')

    # define parameter names for the signal fit
    fit_tpl.SetParName(n_bkgpars, 'N_{sig}')
    fit_tpl.SetParName(n_bkgpars + 1, '#mu')
    fit_tpl.SetParName(n_bkgpars + 2, '#sigma')
    # define parameter values and limits
    fit_tpl.SetParameter(n_bkgpars, 40)
    fit_tpl.SetParLimits(n_bkgpars, -10, 10000)  # modificato
    fit_tpl.SetParameter(n_bkgpars + 1, 2.991)
    fit_tpl.SetParLimits(n_bkgpars + 1, 2.986, 3)

    # define signal and bkg_model TF1 separately
    sigTpl = ROOT.TF1('fitTpl', 'gausn(0)', 0, 5)
    

    # plotting stuff for fit_tpl
    fit_tpl.SetNpx(300)
    fit_tpl.SetLineWidth(2)
    fit_tpl.SetLineColor(2)
    # plotting stuff for bkg model
    bkg_tpl.SetNpx(300)
    bkg_tpl.SetLineWidth(2)
    bkg_tpl.SetLineStyle(2)
    bkg_tpl.SetLineColor(2)

    # define limits for the sigma if provided
    if sigma_limits != None:
        fit_tpl.SetParameter(n_bkgpars + 2, 0.5 *
                             (sigma_limits[0] + sigma_limits[1]))
        fit_tpl.SetParLimits(n_bkgpars + 2, sigma_limits[0], sigma_limits[1])
    # if the mc sigma is provided set the sigma to that value
    elif fixsigma > 0:
        fit_tpl.FixParameter(n_bkgpars + 2, fixsigma)
    # otherwise set sigma limits reasonably
    else:
        fit_tpl.SetParameter(n_bkgpars + 2, 0.002)
        fit_tpl.SetParLimits(n_bkgpars + 2, 0.001, 0.003)

    ########################################
    # plotting the fits
    if mode == 2:
        ax_titles = ';m (^{3}He + #pi) (GeV/#it{c}^{2});Counts' + f' / {round(1000 * histo.GetBinWidth(1), 2)} MeV'
    if mode == 3:
        ax_titles = ';m (d + p + #pi) (GeV/#it{c}^{2});Counts' + f' / {round(1000 * histo.GetBinWidth(1), 2)} MeV'

    # invariant mass distribution histo and fit
    histo.UseCurrentStyle()
    histo.SetLineColor(1)
    histo.SetMarkerStyle(20)
    histo.SetMarkerColor(1)
    histo.SetTitle(ax_titles)
    histo.SetMaximum(1.5 * histo.GetMaximum())
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.SetDrawOption("e")
    histo.GetXaxis().SetRangeUser(2.96, 3.04)
    # represent the bkg_model separately
    bkg_tpl.SetParameters(fit_tpl.GetParameters())
    bkg_tpl.SetLineColor(600)
    bkg_tpl.SetLineStyle(2)
    bkg_tpl.Draw("same")
    # represent the signal model separately
    sigTpl.SetParameter(0, fit_tpl.GetParameter(n_bkgpars))
    sigTpl.SetParameter(1, fit_tpl.GetParameter(n_bkgpars+1))
    sigTpl.SetParameter(2, fit_tpl.GetParameter(n_bkgpars+2))
    sigTpl.SetLineColor(600)
    # sigTpl.Draw("same")

    # get the fit parameters
    mu = fit_tpl.GetParameter(n_bkgpars+1)
    muErr = fit_tpl.GetParError(n_bkgpars+1)
    sigma = fit_tpl.GetParameter(n_bkgpars+2)
    sigmaErr = fit_tpl.GetParError(n_bkgpars+2)
    signal = fit_tpl.GetParameter(n_bkgpars) / histo.GetBinWidth(1)
    errsignal = fit_tpl.GetParError(n_bkgpars) / histo.GetBinWidth(1)
    bkg = bkg_tpl.Integral(mu - nsigma * sigma, mu +
                           nsigma * sigma) / histo.GetBinWidth(1)

    if bkg > 0:
        errbkg = np.sqrt(bkg)
    else:
        errbkg = 0
    # compute the significance
    if signal+bkg > 0:
        signif = signal/np.sqrt(signal+bkg)
        deriv_sig = 1/np.sqrt(signal+bkg)-signif/(2*(signal+bkg))
        deriv_bkg = -signal/(2*(np.power(signal+bkg, 1.5)))
        errsignif = np.sqrt((errsignal*deriv_sig)**2+(errbkg*deriv_bkg)**2)
    else:
        signif = 0
        errsignif = 0

    # print fit info on the canvas
    pinfo2 = ROOT.TPaveText(0.5, 0.5, 0.91, 0.9, "NDC")
    pinfo2.SetBorderSize(0)
    pinfo2.SetFillStyle(0)
    pinfo2.SetTextAlign(30+3)
    pinfo2.SetTextFont(42)

    string = f'ALICE Internal, p-Pb 2018 {cent_class[0]}-{cent_class[1]}%'
    pinfo2.AddText(string)

    decay_label = {
        "": ['{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi^{-} + c.c.', '{}^{3}_{#Lambda}H#rightarrow dp#pi^{-} + c.c.'],
        "_matter": ['{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi^{-}', '{}^{3}_{#Lambda}H#rightarrow dp#pi^{-}'],
        "_antimatter": ['{}^{3}_{#bar{#Lambda}}#bar{H}#rightarrow ^{3}#bar{He}#pi^{+}', '{}^{3}_{#Lambda}H#rightarrow #bar{d}#bar{p}#pi^{+}'],
    }

    string = decay_label[split][mode-2]+', %i #leq #it{ct} < %i cm %i #leq #it{p}_{T} < %i GeV/#it{c} ' % (
        ct_range[0], ct_range[1], pt_range[0], pt_range[1])
    pinfo2.AddText(string)

    string = f'Significance ({nsigma:.0f}#sigma) {signif:.1f} #pm {errsignif:.1f} '
    pinfo2.AddText(string)

    string = f'S ({nsigma:.0f}#sigma) {signal:.0f} #pm {errsignal:.0f}'
    pinfo2.AddText(string)

    string = f'B ({nsigma:.0f}#sigma) {bkg:.0f} #pm {errbkg:.0f}'
    pinfo2.AddText(string)

    if bkg > 0:
        ratio = signal/bkg
        string = f'S/B ({nsigma:.0f}#sigma) {ratio:.4f}'
        pinfo2.AddText(string)

    string = f'#mu {mu:.4f} #pm {muErr:.4f}'
    pinfo2.AddText(string)

    pinfo2.Draw()
    ROOT.gStyle.SetOptStat(0)

    st = histo.FindObject('stats')
    if isinstance(st, ROOT.TPaveStats):
        st.SetX1NDC(0.12)
        st.SetY1NDC(0.62)
        st.SetX2NDC(0.40)
        st.SetY2NDC(0.90)
        st.SetOptStat(0)

    cv.Write()
    
    return (signal, errsignal, signif, errsignif, mu, muErr, sigma, sigmaErr)


def fit_hist_ls(histo, model="pol0", mode=2, split=''):
    # canvas for plotting the invariant mass distribution
    cv = ROOT.TCanvas(f'cv_{histo.GetName()}')

    # define the number of parameters depending on the bkg model
    if 'pol' in str(model):
        n_bkgpars = int(model[3]) + 1
    elif 'expo' in str(model):
        n_bkgpars = 2
    else:
        print(f'Unsupported model {model}')

    # define the fit function bkg_model + gauss
    fit_tpl = ROOT.TF1('fitTpl', f'{model}(0)+gausn({n_bkgpars})', 0, 5)

    # redefine parameter names for the bkg_model
    for i in range(n_bkgpars):
        fit_tpl.SetParName(i, f'B_{i}')

    # define parameter names for the signal fit
    fit_tpl.SetParName(n_bkgpars, 'N_{sig}')
    fit_tpl.SetParName(n_bkgpars + 1, '#mu')
    fit_tpl.SetParName(n_bkgpars + 2, '#sigma')
    # define parameter values and limits
    fit_tpl.SetParameter(n_bkgpars, 40)
    fit_tpl.SetParLimits(n_bkgpars, -10, 10000)  # modificato
    fit_tpl.SetParameter(n_bkgpars + 1, 2.991)
    fit_tpl.SetParLimits(n_bkgpars + 1, 2.986, 3)



    # plotting stuff for fit_tpl
    fit_tpl.SetNpx(300)
    fit_tpl.SetLineWidth(2)
    fit_tpl.SetLineColor(2)
    # plotting stuff for bkg model

    fit_tpl.SetParameter(n_bkgpars + 2, 0.002)
    fit_tpl.SetParLimits(n_bkgpars + 2, 0.001, 0.006)

    ########################################
    # plotting the fits
    if mode == 2:
        ax_titles = ';m (^{3}He + #pi) (GeV/#it{c}^{2});Counts' + f' / {round(1000 * histo.GetBinWidth(1), 2)} MeV'
    if mode == 3:
        ax_titles = ';m (d + p + #pi) (GeV/#it{c}^{2});Counts' + f' / {round(1000 * histo.GetBinWidth(1), 2)} MeV'

    # invariant mass distribution histo and fit
    histo.UseCurrentStyle()
    histo.SetLineColor(1)
    histo.SetMarkerStyle(20)
    histo.SetMarkerColor(1)
    histo.SetTitle(ax_titles)
    histo.SetMaximum(1.5 * histo.GetMaximum())
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.SetDrawOption("e")
    histo.GetXaxis().SetRangeUser(2.96, 3.04)
    # represent the bkg_model separately

    # get the fit parameters
    mu = fit_tpl.GetParameter(n_bkgpars+1)
    muErr = fit_tpl.GetParError(n_bkgpars+1)
    sigma = fit_tpl.GetParameter(n_bkgpars+2)
    sigmaErr = fit_tpl.GetParError(n_bkgpars+2)
    signal = fit_tpl.GetParameter(n_bkgpars) / histo.GetBinWidth(1)
    errsignal = fit_tpl.GetParError(n_bkgpars) / histo.GetBinWidth(1)
    
    cv.Write()

    params = []
    for i in range(n_bkgpars + 3):
        params.append(fit_tpl.GetParameter(i))

    
    return params

def computeAverage(Vals, breakVal = 100):
  Mean = 0
  Stat = 0
  Sys1 = 0
  Sys2 = 0
  weight = 0
  for Bin in Vals:
      w = Bin["bin"][1] - Bin["bin"][0]
      if Bin["bin"][0] >= breakVal:
            break
      weight = weight + w
      Mean = Mean + (Bin["measure"][0] * w)
      Stat = np.hypot(Stat, Bin["measure"][1] * w)
      Sys1 = Sys1 + (Bin["measure"][2] * w)
      Sys2 = Sys2 + (Bin["measure"][3] * w)
  return [Mean / weight, Stat / weight, Sys1 / weight, Sys2 / weight]

def myHypot(x, y, z=0, w=0):
    return np.sqrt(x**2 + y**2 + z**2 + w**2)


def significance_scan(eff_cut, cut_array, dataH, eff_presel, working_point=None, syst_range=0.05):

    hyp_mass = 2.991
    sigma = 0.0015
    m_min = hyp_mass - 3*sigma
    m_max = hyp_mass + 3*sigma

    pol_degree = 0
    bins = 36
    bin_width = (3.04-2.96)/bins
    s_exp = 4.2*85.1*0.25*2*eff_presel

    significance_array = []
    significance_error_array = []

    for i, (bdt_eff, cut) in enumerate(zip(eff_cut, cut_array)):
    
        selected_data_hndl = dataH.get_subset(f"model_output>{cut}")
        counts, bin_edges = np.histogram(selected_data_hndl.get_data_frame()["m"], bins=bins, range=[2.96, 3.04])

        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
        side_map = np.logical_or(bin_centers < m_min, bin_centers > m_max)
        mass_map = np.logical_not(side_map)
        bin_side = bin_centers[side_map]
        bin_mass = bin_centers[mass_map]
        counts_side = counts[side_map]
        counts_mass = counts[mass_map]

        coeff = np.polyfit(bin_side, counts_side, pol_degree)
        data_array = np.arange(2.96, 3.04, bin_width)

        coeff_int = [coeff[0], 0]
        B = (np.polyval(coeff_int, m_max) - np.polyval(coeff_int, m_min))/bin_width
        s = s_exp

        significance_array.append(s/np.sqrt(s+B))
        significance_error_array.append(significance_error(s, B))

    significance = np.asarray(significance_array)
    error_array = np.asarray(significance_error_array)

    significance = significance*eff_cut
    error_array = error_array*eff_cut

    low_limit = significance - error_array
    up_limit = significance + error_array
    fig = plt.figure()
    plt.plot(eff_cut, significance, 'b',   label='Expected Significance')
    plt.fill_between(eff_cut, low_limit, up_limit,
                    facecolor='deepskyblue', label=r'$ \pm 1\sigma$', alpha=0.3)

    if(working_point==None):
        plt.legend(loc="lower center")
    else:
        plt.plot([working_point, working_point],[-1,5.3],linestyle = "--", color = "r", label = "Working Point")
        plt.plot([working_point - syst_range, working_point  - syst_range],[-1,5.3],linestyle = "--", color = "g", label = "Systematic variation range")
        plt.plot([working_point + syst_range, working_point  + syst_range],[-1,5.3],linestyle = "--", color = "g")

        handles, labels = fig.gca().get_legend_handles_labels()
        order = [0,3,1,2]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='lower left')

    plt.xlabel("BDT Efficiency")
    plt.ylabel("Significance x BDT Efficiency")
    plt.xlim(0.5,0.98)
    plt.ylim(0.3, 5.3)
    return fig


def compute_3sigma_coverage(df):
    len_df = len(df)
    sigma_prob = 0.95
    up_val = np.max(df)
    cent_val = np.mean(df)
    var_interval = np.linspace(0, (up_val - cent_val), 1000)
    for var in var_interval:
        len_part = np.sum(np.logical_and(df> cent_val - var, df< cent_val + var))
        frac=  len_part/len_df
        if(frac>sigma_prob):
            return cent_val, var
    