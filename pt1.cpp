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

double gumbel(double* xPtr, double par[]) {
  double x = *xPtr;
  double B = par[0];
  double mu = par[1];
  double A = par[2];
  double z = (x-mu)/B;
  double exponent = -(z + TMath::Exp(-z));
  double f = A*(1/B)*TMath::Exp(exponent);
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
  const int npar2 = 3;

  TMinuit minuit(npar);
  TMinuit minuit2(npar);

  minuit.SetFCN(fcn);
  minuit2.SetFCN(fcn);

  TF1* myfunc = new TF1("myfunc", gausX2, xmin, xmax, npar);
  TF1* mygumbelfunc = new TF1("mygumbelfunc", gumbel, xmin, xmax, npar);

  double par[npar];
  double par2[npar2];
  double stepSize[npar];
  double stepSize2[npar2];
  double minVal[npar];
  double minVal2[npar2];
  double maxVal[npar];
  double maxVal2[npar2];
  TString parName[npar];
  TString parName2[npar2];

  double hmax = h->GetMaximum();
  double hmu = h->GetMean();
  double hstd = h->GetRMS();

  par[0] = 16*hmax;
  par[1] = hmu-2;
  par[2] = hstd-2;
  par[3] = 5*hmax;
  par[4] = hmu+10;
  par[5] = hstd;

  par2[0] = 40;
  par2[1] = 50;
  par2[2] = 1000;

  myfunc->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);

  mygumbelfunc->SetParameters(par[0], par[1], par[2]);
  
  
  
  for (int i=0; i < npar; i++) {
    stepSize[i] = TMath::Abs(par[i]*0.01);
    minVal[i] = par[i]/2;
    maxVal[i] = par[i]*2;
  }

  for (int i=0; i < npar2; i++) {
    stepSize2[i] = TMath::Abs(par2[i]*0.01);
    minVal2[i] = 0;
    maxVal2[i] = 0;
  }

  parName[0] = "A1";
  parName[1] = "mean1";
  parName[2] = "standard deviation 1";
  parName[3] = "A2";
  parName[4] = "mean2";
  parName[5] = "standard deviation 2";

  parName2[0] = "Beta";
  parName2[1] = "mu";
  parName2[2] = "A";

  for (int i=0; i<npar; i++) {
    minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
  }

  cout << "Fail 1" << endl;

  for (int i=0; i<npar2; i++) {
    minuit2.DefineParameter(i, parName2[i].Data(), par2[i], stepSize2[i], minVal2[i], maxVal2[i]);
  }

  cout << "Fail 2" << endl;

  hdata = h;
  fparam = myfunc;

  minuit.Migrad();

  cout << "Fail 3" << endl;

  double outpar[npar];
  double err[npar];
  for (int i=0; i<npar; i++) {
    minuit.GetParameter(i, outpar[i], err[i]);
  }

  cout << "Fail 4" << endl;

  gMinuit->mnmnos();

  myfunc->SetParameters(outpar);

  myfunc->SetLineStyle(1);
  myfunc->SetLineColor(1);
  myfunc->SetLineWidth(2);

  myfunc->GetXaxis()->SetTitle("x");
  myfunc->GetYaxis()->SetTitle("f");

  
  fparam = mygumbelfunc;

  minuit2.Migrad();

  double outpar2[npar2];
  double err2[npar2];
  for (int i=0; i<npar2; i++) {
    minuit2.GetParameter(i, outpar2[i], err2[i]);
  }

  gMinuit->mnmnos();

  mygumbelfunc->SetParameters(outpar2);

  mygumbelfunc->SetLineStyle(1);
  mygumbelfunc->SetLineColor(1);
  mygumbelfunc->SetLineWidth(2);

  mygumbelfunc->GetXaxis()->SetTitle("x");
  mygumbelfunc->GetYaxis()->SetTitle("f");

  

  canvas->Divide(2,1);

  canvas->cd(1);
  h->Draw();
  myfunc->Draw("same");

  canvas->cd(2);
  h->Draw();
  mygumbelfunc->Draw("same");
  
  canvas->SaveAs("pt1.png");

  double chi = 0.0;
  double chi2 = 0.0;
  double measured = 0.0;
  double expected = 0.0;
  double expected2 = 0.0;
  double error = 0.0;
  double binCenter = 0.0;
  for (int i = 0; i < h->GetNbinsX(); i++) {
    measured = h->GetBinContent(i);
    error = h->GetBinError(i);
    binCenter = h->GetBinCenter(i);
    expected = myfunc->Eval(binCenter);
    expected2 = mygumbelfunc->Eval(binCenter);
    if (error>=1e-10) {
      chi += (measured-expected) * (measured-expected) / (error*error);
      chi2 += (measured-expected2) * (measured-expected2) / (error*error);
    }
  }

  cout << "Sum of Two Gaussians:" << endl;
  cout << "chi-squared: " << chi << endl;
  cout << "p value = " << TMath::Prob(chi, h->GetNbinsX()-npar) << endl << endl;

  cout << "Gumbel Distribution:" << endl;
  cout << "chi-squared: " << chi2 << endl;
  cout << "p value = " << TMath::Prob(chi2, h->GetNbinsX()-npar2) << endl;

  cout << "press ENTER to close" << endl;
  cin.get();
 
}
