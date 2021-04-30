using namespace RooStats;
using namespace RooFit;

Double_t Pvalue(Double_t significance) {
  return ROOT::Math::chisquared_cdf_c(pow(significance,2),1)/2;
}

void p0Plot()
{
const int kBlueC  = TColor::GetColor("#1f78b4");
const int kBlueCT = TColor::GetColorTransparent(kBlueC, 0.5);
const int kRedC  = TColor::GetColor("#e31a1c");
const int kRedCT = TColor::GetColorTransparent(kRedC, 0.5);
const int kPurpleC  = TColor::GetColor("#911eb4");
const int kPurpleCT = TColor::GetColorTransparent(kPurpleC, 0.5);
const int kOrangeC  = TColor::GetColor("#ff7f00");
const int kOrangeCT = TColor::GetColorTransparent(kOrangeC, 0.5);
const int kGreenC  = TColor::GetColor("#33a02c");
const int kGreenCT = TColor::GetColorTransparent(kGreenC, 0.5);
const int kMagentaC  = TColor::GetColor("#f032e6");
const int kMagentaCT = TColor::GetColorTransparent(kMagentaC, 0.5);
const int kYellowC  = TColor::GetColor("#ffe119");
const int kYellowCT = TColor::GetColorTransparent(kYellowC, 0.5);
const int kBrownC  = TColor::GetColor("#b15928");

  const char *infile = "../../Utils/Workspaces/ws_eff_0.72.root";
  const char *workspaceName = "ws_eff_0.72";
  const char *modelConfigName = "ModelConfig";
  const char *dataName = "data";

  /////////////////////////////////////////////////////////////
  // First part is just to access the workspace file
  ////////////////////////////////////////////////////////////

  // Check if example input file exists
  TFile *file = TFile::Open(infile);

  // get the workspace out of the file
  RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);

  // get the modelConfig out of the file
  RooStats::ModelConfig *mc = (RooStats::ModelConfig *)w->obj(modelConfigName);

  // get the modelConfig out of the file
  RooAbsData *data = w->data(dataName);

  // get the modelConfig (S+B) out of the file
  // and create the B model from the S+B model
  ModelConfig *sbModel = (RooStats::ModelConfig *)w->obj(modelConfigName);
  sbModel->SetName("S+B Model");
  RooRealVar *poi = (RooRealVar *)sbModel->GetParametersOfInterest()->first();
  poi->Print();
  sbModel->SetSnapshot(*poi);
  ModelConfig *bModel = (ModelConfig *)sbModel->Clone();
  bModel->SetName("B Model");
  poi->setVal(0.);
  bModel->SetSnapshot(*poi);
  bModel->Print();

  vector<double> masses;
  vector<double> p0values;
  vector<double> p0valuesExpected;

  double massMin = -0.06;
  double massMax = 0.06;

  // loop on the mass values
  

  for (double mass = massMin; mass <= massMax; mass += (massMax - massMin) / 400)
  {
    cout << endl
         << endl
         << "Running for mass: " << 2.9901052 + mass << endl
         << endl;
    w->var("deltaM")->setVal(mass);
    w->var("deltaM")->setConstant(true);
    w->var("deltaM")->Print();

    AsymptoticCalculator *ac = new AsymptoticCalculator(*data, *sbModel, *bModel);
    ac->SetOneSidedDiscovery(true); // for one-side discovery test
    AsymptoticCalculator::SetPrintLevel(-1);

    HypoTestResult *asymCalcResult = ac->GetHypoTest();

    asymCalcResult->Print();
    w->var("deltaM")->Print();
    masses.push_back(2.9901052 - mass);
    p0values.push_back(asymCalcResult->NullPValue());

    double expectedP0 = AsymptoticCalculator::GetExpectedPValues(asymCalcResult->NullPValue(), asymCalcResult->AlternatePValue(), 0, false);
    p0valuesExpected.push_back(expectedP0);
    std::cout << "expected p0 = " << expectedP0 << std::endl;
  }


  TFile fit_file("../../Results/inv_mass_fits.root");
  RooPlot* frame = static_cast<RooPlot*>(fit_file.Get("frame_0.72"));
  

  constexpr double kYsize = 0.4;

  gStyle->SetOptStat(1);
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetLabelSize(0.04,"xyz"); // size of axis value font
  gStyle->SetTitleSize(0.05,"xyz"); // size of axis title font
  gStyle->SetTitleFont(42,"xyz"); // font option
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(1.1,"y");
  // default canvas options
  gStyle->SetCanvasDefW(800);
  gStyle->SetCanvasDefH(600);
  gStyle->SetPadBottomMargin(0.12); //margins...
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadGridX(0); // grids, tickmarks
  gStyle->SetPadGridY(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPaperSize(20,24); // US letter size
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);


  TPaveText* pinfo_alice = new TPaveText(0.45, 0.75, 0.93, 0.86, "NDC");
  pinfo_alice->SetBorderSize(0);
  pinfo_alice->SetFillStyle(0);
  pinfo_alice->SetTextFont(42);
  pinfo_alice->AddText("ALICE");
  pinfo_alice->AddText("p-Pb 0-40%, #sqrt{#it{s}_{NN}} = 5.02 TeV");

  TCanvas *cv = new TCanvas("cv1","cv1",1500,1500);
  cv->cd();

  
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitle("Entries / (2.35 MeV/#it{c}^{2})");
  frame->GetXaxis()->SetTitle("{}_{#Lambda}^{3}H mass (^{3}He + #pi) (GeV/#it{c}^{2})");
  frame->SetMinimum(0.01);
  frame->Draw();
  pinfo_alice->Draw();
  frame->Print();
  TLegend *leg1 = new TLegend(0.49,0.56,0.9,0.75);
  leg1->AddEntry("h_data","{}_{#Lambda}^{3}H + {}_{#bar{#Lambda}}^{3}#bar{H}", "PE");
  leg1->AddEntry("model_Norm[m]_Range[fit_nll_model_data]_NormRange[fit_nll_model_data]","Signal + Background", "L");

  leg1->AddEntry("model_Norm[m]_Comp[bkg]_Range[fit_nll_model_data]_NormRange[fit_nll_model_data]","Background", "L");
  leg1->Draw();

  cv->SaveAs("../../Results/signal_extraction.png");
  cv->SaveAs("../../Results/signal_extraction.pdf");




  TGraph *graph1 = new TGraph(masses.size(), &masses[0], &p0values[0]);
  TGraph *graph2 = new TGraph(masses.size(), &masses[0], &p0valuesExpected[0]);

  TLine *l1 = new TLine(2.96, Pvalue(1), 3.04, Pvalue(1));
  TLine *l2 = new TLine(2.96, Pvalue(2), 3.04, Pvalue(2));
  TLine *l3 = new TLine(2.96, Pvalue(3), 3.04, Pvalue(3));
  TLine *l4 = new TLine(2.96, Pvalue(4), 3.04, Pvalue(4));
  TLine *l5 = new TLine(2.96, Pvalue(5), 3.04, Pvalue(5));


  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l3->SetLineStyle(2);
  l4->SetLineStyle(2);
  l5->SetLineStyle(2);

  TPaveText* pinfo2 = new TPaveText(0.85, 0.73, 0.91, 0.82, "NDC");
  pinfo2->SetBorderSize(0);
  pinfo2->SetFillStyle(0);
  pinfo2->SetTextAlign(30+3);
  pinfo2->SetTextFont(42);
  pinfo2->AddText("2#sigma");

  TPaveText* pinfo3 = new TPaveText(0.85, 0.6, 0.91, 0.69, "NDC");
  pinfo3->SetBorderSize(0);
  pinfo3->SetFillStyle(0);
  pinfo3->SetTextAlign(30+3);
  pinfo3->SetTextFont(42);
  pinfo3->AddText("3#sigma");

  TPaveText* pinfo4 = new TPaveText(0.85, 0.425,0.91, 0.515, "NDC");
  pinfo4->SetBorderSize(0);
  pinfo4->SetFillStyle(0);
  pinfo4->SetTextAlign(30+3);
  pinfo4->SetTextFont(42);
  pinfo4->AddText("4#sigma");

  TPaveText* pinfo5 = new TPaveText(0.85, 0.21,0.91, 0.3, "NDC");
  pinfo5->SetBorderSize(0);
  pinfo5->SetFillStyle(0);
  pinfo5->SetTextAlign(30+3);
  pinfo5->SetTextFont(42);
  pinfo5->AddText("5#sigma");




  graph1->SetMarkerStyle(20);


  TCanvas *cv2 = new TCanvas("cv2","cv2",1500,1500);
  cv2->cd();



  cv2->SetLogy();
  graph1->GetYaxis()->SetTitleSize(0.06);
  graph1->GetYaxis()->SetTitleOffset(0.9);
  graph1->GetXaxis()->SetTitleOffset(1.1);
  graph1->GetXaxis()->SetTitle("{}_{#Lambda}^{3}H mass (^{3}He + #pi) (GeV/#it{c}^{2})");
  graph1->GetYaxis()->SetTitle("Local p_{0}");
  graph1->GetXaxis()->SetLimits(2.96,3.04);
  graph1->GetXaxis()->SetRangeUser(2.96,3.04);
  graph1->Draw("AC same");
  
  graph1->SetTitle("");
  graph1->SetMaximum(0.9);
  graph1->SetMinimum(0.5e-7);
  graph1->SetLineColor(kBlueCT);
  graph1->SetLineWidth(4);


  l2->Draw();
  l3->Draw();
  l4->Draw();
  l5->Draw();
  pinfo2->Draw();
  pinfo3->Draw();
  pinfo4->Draw();
  pinfo5->Draw();

  cv2->SaveAs("../../Results/significance.png");
  cv2->SaveAs("../../Results/significance.pdf");

  TFile sig_file("signal_extraction.root", "recreate");
  cv->Write();
  cv2->Write();
  sig_file.Close();

}
