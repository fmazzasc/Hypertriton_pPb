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

void GenerateTableFromData(bool likeSign = false, bool pp = false)
{

  string dataDir = "../Trees/";
  string tableDir = "../Tables/";

  string lsString = likeSign ? "_LS.root" : ".root";

  string inFileName16 = pp ? "HyperTritonTree_pp13TeV" : "HyperTritonTree_16qt";

  
  string inFileArg16 = dataDir + inFileName16 + lsString;

  string inFileName13 = "HyperTritonTree_13bc";
  string inFileArg13 = dataDir  + inFileName13 + lsString;

  string outFileName = pp ? "DataTable_pp" : "DataTable_pPb";
  string outFileArg = tableDir  + outFileName + lsString;

  string treeName = pp ? "Hyp3O2" : "_default/fTreeV0";

  TChain inputChain(treeName.data());
  inputChain.AddFile(inFileArg16.data());
  if (!likeSign && !pp)
    inputChain.AddFile(inFileArg13.data());

  TTreeReader fReader(&inputChain);
  TFile outFile(outFileArg.data(), "RECREATE");

  if(!pp){
  TTreeReaderArray<RHyperTritonHe3pi> RHyperVec = {fReader, "RHyperTriton"};
  TTreeReaderValue<RCollision> RColl = {fReader, "RCollision"}; 
  Table2 tree("DataTable", "Data Table");
  while (fReader.Next())
  {
    for (auto &RHyper : RHyperVec)
      tree.Fill(RHyper, *RColl);
  }
  outFile.cd();
  tree.Write();
  }

  if(pp){
  TTreeReaderValue<RHyperTriton3O2> RHyper = {fReader, "RHyperTriton"};
  TableO2 tree(false);
    while (fReader.Next()) {
      tree.Fill(*RHyper);
  }
  outFile.cd();
  tree.Write();
  }


  outFile.Close();
  std::cout << "\nDerived tree from Data generated!\n" << std::endl;

}
