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

TH1F *hdata1;
TH1F *hdata2;
TF1 *fparam1;
TF1 *fparam2;

double func1(double* xPtr, double par[]) {
  double x = *xPtr;
  double A = par[0];
  double mu = par[1];
  double std = par[2];
  double gaus = A*TMath::Exp(-.5*((x-mu)/std)*((x-mu)/std))/(std*TMath::Sqrt(2*TMath::Pi()));
  double lam = par[3];
  double B = par[4];
  double exponential = B * TMath::Exp(-x/lam);
  double f = gaus + exponential;
  return f;
}

double func2(double* xPtr, double par[]) {
  double x = *xPtr;
  double A = par[0];
  double mu = par[1];
  double std = par[2];
  double gaus = A*TMath::Exp(-.5*((x-mu)/std)*((x-mu)/std))/(std*TMath::Sqrt(2*TMath::Pi()));
  double n = par[3];
  double B = par[4];
  double exponential = B * TMath::Power(x, n);
  double f = gaus + exponential;
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

  for (int i=0; i<npar; i++){ //both functions share a gausian
    if (i<3) {
      fparam1->SetParameter(i,par[i]);
      fparam2->SetParameter(i,par[i]);
    } else if (i>=3 && i<5) { //these params are for the falling exponential
      fparam1->SetParameter(i,par[i]);
    } else {
      fparam2->SetParameter(i-2, par[i]);
    }
  }

  f = calcCHI(hdata1,fparam1) + calcCHI(hdata2,fparam2);
 
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

  TFile* myfile = new TFile("experiments.root");
  TH1F* h1 = (TH1F*) myfile->Get("hexp1");
  TH1F* h2 = (TH1F*) myfile->Get("hexp2");
  
  
  
  double xmin = 0.0;
  double xmax = 100.0;

  const int npar1 = 5;
  const int npar2 = 5;
  const int npar = npar1 + npar2 - 3;

  TMinuit minuit(npar);

  minuit.SetFCN(fcn);

  TF1* myfunc1 = new TF1("myfunc1", func1, xmin, xmax, npar1);
  TF1* myfunc2 = new TF1("myfunc2", func2, xmin, xmax, npar2);

  double par[npar];
  double stepSize[npar];
  double minVal[npar];
  double maxVal[npar];
  TString parName[npar];

  par[0] = 500;
  par[1] = 75;
  par[2] = 4.5;
  par[3] = 20;
  par[4] = 700;
  par[5] = -2.2;
  par[6] = 500000;

  myfunc1->SetParameters(par[0], par[1], par[2], par[3], par[4]);
  myfunc2->SetParameters(par[0], par[1], par[2], par[5], par[6]);
  
  
  
  for (int i=0; i < npar; i++) {
    stepSize[i] = TMath::Abs(par[i]*0.01);
    minVal[i] = 0;
    maxVal[i] = 0;
  }

  parName[0] = "A";
  parName[1] = "mean";
  parName[2] = "standard deviation";
  parName[3] = "lambda";
  parName[4] = "func1 multiplier";
  parName[5] = "n";
  parName[6] = "func2 multiplier";

  for (int i=0; i<npar; i++) {
    minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
  }

  hdata1 = h1;
  hdata2 = h2;
  fparam1 = myfunc1;
  fparam2 = myfunc2;

  minuit.Migrad();

  double outpar[npar];
  double err[npar];
  for (int i=0; i<npar; i++) {
    minuit.GetParameter(i, outpar[i], err[i]);
  }

  gMinuit->mnmnos();

  myfunc1->SetParameters(outpar[0], outpar[1], outpar[2], outpar[3], outpar[4]);
  myfunc2->SetParameters(outpar[0], outpar[1], outpar[2], outpar[5], outpar[6]);


  double chi = 0.0;

  double measured1 = 0.0;
  double expected1 = 0.0;
  double error1 = 0.0;
  double binCenter1 = 0.0;
  for (int i = 0; i < h1->GetNbinsX(); i++) {
    measured1 = h1->GetBinContent(i);
    error1 = h1->GetBinError(i);
    binCenter1 = h1->GetBinCenter(i);
    expected1 = myfunc1->Eval(binCenter1);
    if (error1>=1e-10) {
      chi += (measured1-expected1) * (measured1-expected1) / (error1*error1);
    }
  }

  double measured2 = 0.0;
  double expected2 = 0.0;
  double error2 = 0.0;
  double binCenter2 = 0.0;
  for (int i = 0; i < h2->GetNbinsX(); i++) {
    measured2 = h2->GetBinContent(i);
    error2 = h2->GetBinError(i);
    binCenter2 = h2->GetBinCenter(i);
    expected2 = myfunc2->Eval(binCenter2);
    if (error2>=1e-10) {
      chi += (measured2-expected2) * (measured2-expected2) / (error2*error2);
    }
  }


  cout << "chi-squared: " << chi << endl;
  cout << "p value = " << TMath::Prob(chi, h1->GetNbinsX() + h2->GetNbinsX() - npar) << endl;
  

  canvas->Divide(2,1);
  
  canvas->cd(1);
  h1->Draw();
  myfunc1->Draw("Same");
  canvas->cd(2);
  h2->Draw();
  myfunc2->Draw("Same");

  canvas->SaveAs("pt2.png");

  cout << "press ENTER to close" << endl;
  cin.get();
 
}
