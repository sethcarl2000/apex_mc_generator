#include "TROOT.h"
#include "TH1D.h" 

extern "C" void my_subrt_(float *x,float *y,float *z,float *ans);

using namespace std; 

void make_test_histogram() {

  auto hist = new TH1D("h", "Test histogram", 200, -3, 3);

  for (int i=0; i<2e6; i++) {

    float x,y,z,ans;

    my_subrt_(&x, &y, &z, &ans);

    hist->Fill(ans);

  }

  hist->Draw();

}

 
