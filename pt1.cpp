#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <TMinuit.h>
#include <TRandom2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>

using namespace std;

TH1F *hdata;
TF1 *fparam;

double gausX2(double* xPtr, double par[]) {
  double x = *xPtr;
  double A1 = par[0];
  double mu1 = par[1];
  double std1 = par[2];
  double gaus1 = A1*TMath::Exp(-.5*((x-mu1)/std1)*((x-mu1)/std1))/(std1*TMath::Sqrt(2*TMath::Pi()));
  double A2 = par[3];
  double mu2 = par[4];
  double std2 = par[5];
  double gaus2 = A2*TMath::Exp(-.5*((x-mu2)/std2)*((x-mu2)/std2))/(std2*TMath::Sqrt(2*TMath::Pi()));
  double f = gaus1+gaus2;
  return f;
}

double calcCHI(TH1F* h, TF1* f){
  double chi=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    double x=h->GetBinCenter(i);
    double n=h->GetBinContent(i);
    double mu=f->Eval(x);
    double err = h->GetBinError(i);
    if (err>=1e-10){
      chi += (n-mu) * (n-mu) / (err*err);
    }
  }
  return chi;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  for (int i=0; i<npar; i++){
    fparam->SetParameter(i,par[i]);
  }

  f = calcCHI(hdata,fparam);
 
}

int main(int argc, char **argv) {

  TApplication theApp("App", &argc, argv);

  TCanvas* canvas = new TCanvas();

  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();

  TFile* myfile = new TFile("distros.root");
  TH1F* h = (TH1F*) myfile->Get("dist1");
  
  double xmin = 0.0;
  double xmax = 200;

  const int npar = 6;

  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  TF1* myfunc = new TF1("myfunc", gausX2, xmin, xmax, npar);

  double par[npar];
  double stepSize[npar];
  double minVal[npar];
  double maxVal[npar];
  TString parName[npar];

  double hmax = h->GetMaximum();
  double hmu = h->GetMean();
  double hstd = h->GetRMS();

  par[0] = 16*hmax;
  par[1] = hmu-2;
  par[2] = hstd-2;
  par[3] = 5*hmax;
  par[4] = hmu+10;
  par[5] = hstd;

  myfunc->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);
  
  
  
  for (int i=0; i < npar; i++) {
    stepSize[i] = TMath::Abs(par[i]*0.01);
    minVal[i] = par[i]/2;
    maxVal[i] = par[i]*2;
  }

  parName[0] = "A1";
  parName[1] = "mean1";
  parName[2] = "standard deviation 1";
  parName[3] = "A2";
  parName[4] = "mean2";
  parName[5] = "standard deviation 2";

  for (int i=0; i<npar; i++) {
    minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
  }

  hdata = h;
  fparam = myfunc;

  minuit.Migrad();

  double outpar[npar];
  double err[npar];
  for (int i=0; i<npar; i++) {
    minuit.GetParameter(i, outpar[i], err[i]);
  }

  gMinuit->mnmnos();

  myfunc->SetParameters(outpar);

  myfunc->SetLineStyle(1);
  myfunc->SetLineColor(1);
  myfunc->SetLineWidth(2);

  myfunc->GetXaxis()->SetTitle("x");
  myfunc->GetYaxis()->SetTitle("f");


  h->Draw();
  myfunc->Draw("same");
  canvas->SaveAs("pt1.png");

  cout << "press ENTER to close" << endl;
  cin.get();
 
}
