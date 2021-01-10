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

  const char *infile = "../../Utils/Workspaces/ws_pol1_72.root";
  const char *workspaceName = "ws_pol1_72";
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

  double massMin = 2.96;
  double massMax = 3.04;

  // loop on the mass values

  for (double mass = massMin; mass <= massMax; mass += (massMax - massMin) / 40)
  {
    cout << endl
         << endl
         << "Running for mass: " << mass << endl
         << endl;
    w->var("hyp_mass")->setVal(mass);
    w->var("hyp_mass")->setConstant(true);
    w->var("hyp_mass")->Print();

    AsymptoticCalculator *ac = new AsymptoticCalculator(*data, *sbModel, *bModel);
    ac->SetOneSidedDiscovery(true); // for one-side discovery test
    AsymptoticCalculator::SetPrintLevel(-1);

    HypoTestResult *asymCalcResult = ac->GetHypoTest();

    asymCalcResult->Print();
    w->var("hyp_mass")->Print();
    masses.push_back(mass);
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

  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  l3->SetLineStyle(2);
  l4->SetLineStyle(2);


  TPaveText* pinfo2 = new TPaveText(0.94, 0.77, 1, 0.86, "NDC");
  pinfo2->SetBorderSize(0);
  pinfo2->SetFillStyle(0);
  pinfo2->SetTextAlign(30+3);
  pinfo2->SetTextFont(42);
  pinfo2->AddText("2#sigma");

  TPaveText* pinfo3 = new TPaveText(0.94, 0.61, 1, 0.7, "NDC");
  pinfo3->SetBorderSize(0);
  pinfo3->SetFillStyle(0);
  pinfo3->SetTextAlign(30+3);
  pinfo3->SetTextFont(42);
  pinfo3->AddText("3#sigma");

  TPaveText* pinfo4 = new TPaveText(0.94, 0.39,1, 0.48, "NDC");
  pinfo4->SetBorderSize(0);
  pinfo4->SetFillStyle(0);
  pinfo4->SetTextAlign(30+3);
  pinfo4->SetTextFont(42);
  pinfo4->AddText("4#sigma");

  graph1->SetMarkerStyle(20);

  TFile fit_file("../../Results/inv_mass_fits.root");
  TCanvas* canv = static_cast<TCanvas*>(fit_file.Get("cv_templ_0.72_pol1_040"));

  TCanvas *cv = new TCanvas("cv","cv");
  cv->Divide(1,2,0,0);
  auto pad1 = cv->cd(1);
  canv->DrawClonePad();
  pad1->SetBottomMargin(0);
  canv->SetTopMargin(0.);
  pad1->SetTopMargin(0.);

  auto pad2 = cv->cd(2);
  pad2->SetLogy();
  pad2->SetRightMargin(0.05);
  pad2->SetLeftMargin(0.12);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.18);
  graph1->Draw("AC");
  graph1->GetXaxis()->SetTitle("{}_{#Lambda}^{3}H mass (GeV/#it{c}^{2})");
  graph1->GetYaxis()->SetTitle("Local p_{0}");
  graph1->GetXaxis()->SetLimits(2.96,3.04);
  graph1->GetXaxis()->SetRangeUser(2.96,3.04);
  
  graph1->SetTitle("");
  graph1->SetMaximum(0.9);
  graph1->SetMinimum(5e-7);
  graph1->SetLineColor(kBlueCT);
  graph1->SetLineWidth(3);
  graph1->GetXaxis()->SetLabelOffset(0.01);
  graph1->GetXaxis()->SetNdivisions(510);
  graph1->GetXaxis()->SetTitleSize(0.07);
  graph1->GetYaxis()->SetTitleOffset(0.5);
  graph1->GetYaxis()->SetTitleSize(0.07);

  l2->Draw();
  l3->Draw();
  l4->Draw();
  pinfo2->Draw();
  pinfo3->Draw();
  pinfo4->Draw();


  cv->SaveAs("../../Results/significance.png");
  cv->SaveAs("../../Results/significance.pdf");
  cv->SaveAs("../../Results/significance.root");
}


