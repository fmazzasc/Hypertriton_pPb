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

  TPaveText* pinfo2 = new TPaveText(0.94, 0.81, 1, 0.9, "NDC");
  pinfo2->SetBorderSize(0);
  pinfo2->SetFillStyle(0);
  pinfo2->SetTextAlign(30+3);
  pinfo2->SetTextFont(42);
  pinfo2->AddText("2#sigma");

  TPaveText* pinfo3 = new TPaveText(0.94, 0.67, 1, 0.76, "NDC");
  pinfo3->SetBorderSize(0);
  pinfo3->SetFillStyle(0);
  pinfo3->SetTextAlign(30+3);
  pinfo3->SetTextFont(42);
  pinfo3->AddText("3#sigma");

  TPaveText* pinfo4 = new TPaveText(0.94, 0.49,1, 0.58, "NDC");
  pinfo4->SetBorderSize(0);
  pinfo4->SetFillStyle(0);
  pinfo4->SetTextAlign(30+3);
  pinfo4->SetTextFont(42);
  pinfo4->AddText("4#sigma");

  TPaveText* pinfo5 = new TPaveText(0.94, 0.27,1, 0.36, "NDC");
  pinfo5->SetBorderSize(0);
  pinfo5->SetFillStyle(0);
  pinfo5->SetTextAlign(30+3);
  pinfo5->SetTextFont(42);
  pinfo5->AddText("5#sigma");


  TPaveText* pinfo_alice = new TPaveText(0.5, 0.574, 0.92, 0.7, "NDC");
  pinfo_alice->SetBorderSize(0);
  pinfo_alice->SetFillStyle(0);
  pinfo_alice->SetTextFont(42);
  pinfo_alice->AddText("ALICE Internal");
  pinfo_alice->AddText("p-Pb 0-40%, #sqrt{#it{s}_{NN}}=5.02 TeV");

  graph1->SetMarkerStyle(20);

  TFile fit_file("../../Results/inv_mass_fits.root");
  RooPlot* frame = static_cast<RooPlot*>(fit_file.Get("frame_0.72"));
  

  constexpr double kYsize = 0.4;

  TCanvas *cv = new TCanvas("cv","cv",1200,1200);
  TPad *pad0 = new TPad("pad0","",0, kYsize, 1, 1);
  // auto pads = CreatePads(cv);
  // pads[0]->cd();
  cv->cd();
  pad0->SetBottomMargin(0.02);
  pad0->Draw();
  pad0->cd();
  
  frame->GetYaxis()->SetTickLength(0.012 / (1. - kYsize));
  frame->GetYaxis()->SetTitleSize(40);
  frame->GetYaxis()->SetTitleFont(43);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  frame->GetYaxis()->SetLabelSize(24);
  frame->GetXaxis()->SetLabelOffset(100);
  frame->GetXaxis()->SetTitleOffset(100);
  
  frame->Draw();
  pinfo_alice->Draw();
  

  TPad *pad1 = new TPad("pad1","",0, 0, 1, kYsize);
  cv->cd();
  pad1->SetTopMargin(0);
  pad1->SetBottomMargin(0.2);
  pad1->Draw();
  pad1->cd();


  pad1->SetLogy();

  graph1->GetYaxis()->SetTickLength(0.012 / kYsize);
  graph1->GetYaxis()->SetTitleSize(40);
  graph1->GetYaxis()->SetTitleFont(43);
  graph1->GetYaxis()->SetTitleOffset(1.4);
  graph1->GetYaxis()->SetLabelOffset(0.01);
  graph1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  graph1->GetYaxis()->SetLabelSize(24);
  graph1->GetXaxis()->SetTickLength(0.02 / kYsize);
  graph1->GetXaxis()->SetTitleSize(40);
  graph1->GetXaxis()->SetTitleFont(43);
  graph1->GetXaxis()->SetTitleOffset(1.05 / kYsize);
  graph1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  graph1->GetXaxis()->SetLabelSize(24);
  graph1->GetXaxis()->SetTitle("{}_{#Lambda}^{3}H mass (^{3}He + #pi) (GeV/#it{c}^{2})");
  graph1->GetYaxis()->SetTitle("Local p_{0}");
  graph1->GetXaxis()->SetLimits(2.96,3.04);
  graph1->GetXaxis()->SetRangeUser(2.96,3.04);
  graph1->Draw("AC same");
  
  graph1->SetTitle("");
  graph1->SetMaximum(0.9);
  graph1->SetMinimum(0.5e-7);
  graph1->SetLineColor(kBlueCT);
  graph1->SetLineWidth(3);


  l2->Draw();
  l3->Draw();
  l4->Draw();
  l5->Draw();
  pinfo2->Draw();
  pinfo3->Draw();
  pinfo4->Draw();
  pinfo5->Draw();


  cv->SaveAs("../../Results/significance.png");
  cv->SaveAs("../../Results/significance.pdf");
  cv->SaveAs("../../Results/significance.root");
}
