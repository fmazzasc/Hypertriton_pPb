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

void GenerateTableFromData(bool likeSign = false)
{

  string dataDir = "../Trees/";
  string tableDir = "../Tables/";

  string lsString = likeSign ? "_LS.root" : ".root";

  string inFileName16 = "HyperTritonTree_16qt";
  string inFileArg16 = dataDir + inFileName16 + lsString;

  string inFileName13 = "HyperTritonTree_13bc";
  string inFileArg13 = dataDir  + inFileName13 + lsString;

  string outFileName = "DataTable_pPb";
  string outFileArg = tableDir  + outFileName + lsString;

  TChain inputChain("_default/fTreeV0");
  inputChain.AddFile(inFileArg16.data());
  if (!likeSign)
    inputChain.AddFile(inFileArg13.data());

  TTreeReader fReader(&inputChain);
  TTreeReaderArray<RHyperTritonHe3pi> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"};

  TFile outFile(outFileArg.data(), "RECREATE");
  Table2 tree("DataTable", "Data Table");

  while (fReader.Next())
  {
    for (auto &RHyper : RHyperVec)
      tree.Fill(RHyper, *RColl);
  }

  outFile.cd();
  tree.Write();
  outFile.Close();

  std::cout << "\nDerived tree from Data generated!\n"
            << std::endl;
}
