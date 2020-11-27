#include <cmath>
#include <string>
#include "AliExternalTrackParam.h"
#include "ROOT/RDataFrame.hxx"
#include "TDirectory.h"
#include "TF1.h"
#include "TH1D.h"

double BetheBlochAleph(double *x, double *p) {
  return AliExternalTrackParam::BetheBlochAleph(x[0] / 2.80923f, p[0], p[1], p[2], p[3], p[4]);
}

void CalibHe3(TH1D* mean, TH1D *sigma, TFile* tfile) {

  sigma->Divide(mean);
  sigma->Fit("pol0");

  TF1 mybethe("mybethe",BetheBlochAleph, 1, 6, 5);
  double starting_pars[5]{-166.11733,-0.11020473,0.10851357,2.7018593,-0.087597824};
  mybethe.SetParameters(starting_pars);

  for (int i{0}; i < 10; ++i)
    mean->Fit(&mybethe);
  tfile->cd();
  mean->Write();
  sigma->Write();

  std::cout<< "par0: " << mybethe.GetParameter(0) << std::endl;
  std::cout<< "par1: " << mybethe.GetParameter(1) << std::endl;
  std::cout<< "par2: " << mybethe.GetParameter(2) << std::endl;
  std::cout<< "par3: " << mybethe.GetParameter(3) << std::endl;
  std::cout<< "par4: " << mybethe.GetParameter(4) << std::endl;

}

