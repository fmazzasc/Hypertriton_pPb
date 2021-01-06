#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "AliAnalysisTaskHyperTriton2He3piML.h"
#include "AliPID.h"

#include "helpers/Common.h"
#include "helpers//GenTable2.h"
#include "helpers/Table2.h"

void GenerateTableFromMC(bool reject = true, string ptShape = "mtexp")
{
  gRandom->SetSeed(1995);

  string inFileName = "HyperTritonTree_17d.root";
  string inFileArg =  "../Trees/" + inFileName;

  string outFileName = "SignalTable_17d_" + ptShape + ".root";
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
  

  float max = hypPtShape->GetMaximum();

  TFile *inFile = new TFile(inFileArg.data(), "READ");
  TTreeReader fReader("_default/fTreeV0", inFile);
  TTreeReaderArray<RHyperTritonHe3pi> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderArray<SHyperTritonHe3pi> SHyperVec = {fReader, "SHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"};


  // new flat tree with the features
  TFile outFile(outFileArg.data(), "RECREATE");
  Table2 table("SignalTable", "Signal Table");
  GenTable2 genTable("GenTable", "Generated particle table");

  while (fReader.Next())
  {
    auto cent = RColl->fCent;
    for (auto &SHyper : SHyperVec)
    {

      bool matter = SHyper.fPdgCode > 0;

      double pt = std::hypot(SHyper.fPxHe3 + SHyper.fPxPi, SHyper.fPyHe3 + SHyper.fPyPi);
      
      if (reject)
      {
        float hypPtShapeNum = hypPtShape->Eval(pt) / max;
        if (hypPtShapeNum < gRandom->Rndm())
          continue;
      }
      genTable.Fill(SHyper, *RColl);
      int ind = SHyper.fRecoIndex;

      if (ind >= 0)
      {
        auto &RHyper = RHyperVec[ind];
        table.Fill(RHyper, *RColl);
        double recpt = std::hypot(RHyper.fPxHe3 + RHyper.fPxPi, RHyper.fPyHe3 + RHyper.fPyPi);
        
      }
    }
  }

  outFile.cd();

  table.Write();
  genTable.Write();
  outFile.Close();


}
