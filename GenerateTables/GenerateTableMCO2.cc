#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

using namespace std;

#include "helpers/Common.h"
#include "helpers/Table2.h"
#include "helpers/O2Table.h"


void GenerateTableMCO2(bool reject = true, string ptShape = "mtexp")
{

  gRandom->SetSeed(1995);
  string inFileName = "HyperTritonTree_pp13TeV_MC.root";
  string inFileArg =  "../Trees/" + inFileName;
  string outFileName = "SignalTable_pp13TeV_" + ptShape + ".root";

  string bwFileName;
  string outFileArg =  "../Tables/" + outFileName;

  if(ptShape=="bw")
    bwFileName = "output_He3_yieldfits.root";
  else
    bwFileName = "FittedFunctions_He3.root";

  string bwFileArg = "helpers/" + bwFileName;

  TFile bwFile(bwFileArg.data());

  TF1 *hypPtShape{nullptr};

  if(ptShape=="mtexp")
    hypPtShape = (TF1 *)bwFile.Get("#it{m}_{T} exponentialfittedtoMultiplicityIntegrated");
  else if(ptShape=="ptexp")
    hypPtShape = (TF1 *)bwFile.Get("#it{p}_{T} exponentialfittedtoMultiplicityIntegrated");
  else if(ptShape=="bol")
    hypPtShape = (TF1 *)bwFile.Get("BoltzmannfittedtoMultiplicityIntegrated");
  else if(ptShape=="bw")
    hypPtShape = (TF1 *)bwFile.Get("He3/BlastWave/BlastWave0");

  else
  {
    std::cout << ptShape << " not found" << std::endl;
    return;
  }

  TChain inputChain("Hyp3O2");
  inputChain.Add(inFileArg.data());

  TFile outFile(outFileArg.data(), "RECREATE");
  TableO2 tree(true);

  TTreeReader fReader(&inputChain);

TTreeReaderValue<SHyperTriton3O2> SHyper{fReader, "SHyperTriton"};

while (fReader.Next()) {
    tree.Fill(*SHyper, hypPtShape, hypPtShape->GetMaximum());
}

outFile.cd();
tree.Write();
outFile.Close();

std::cout << "\nDerived tables from MC generated!\n" << std::endl;

}