#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include <cmath>

int plotRmsFits()
 {
  auto rd=TFile::Open("triple-run068-lgad04-07-08-3GeV-40pct-v1.root","READ");
  if (!rd->IsOpen())
   return 1;
  auto mean=new TH1F("baselineRMSmean","Baseline RMS mean value;Channel (number);RMS fitted mean (mV)",30,-0.5,29.5);
  auto sig=new TH1F("baselineRMSsig","Baseline RMS fitted sigma;Channel (number);RMS fitted sigma (mV)",30,-0.5,29.5);
  auto sixsig=new TH1F("baselineRMS6sig","Baseline RMS 6-sigma cut;Channel (number);RMS 6 sigma cut (mV)",30,-0.5,29.5);
  TCanvas c1;
  gStyle->SetOptFit(1);
  for (uint8_t ch=0; ch < 30; ch++)
   {
    std::string hname="histBaselineVar_ch"+std::to_string(ch)+";1";
    c1.Clear();
    auto rms=rd->Get<TH1F>(hname.c_str());
    auto resp=rms->Fit("gaus","S","",0.,2.);
    mean->SetBinContent(ch+1,resp->Parameter(1));
    mean->SetBinError(ch+1,resp->ParError(1));
    sig->SetBinContent(ch+1,resp->Parameter(2));
    sig->SetBinError(ch+1,resp->ParError(2));
    auto sig6=resp->Parameter(1)+6.*resp->Parameter(2);
    auto s6m=resp->ParError(1);
    auto s6s=6.*resp->ParError(2);
    auto s6e=std::sqrt(s6s*s6s+s6m*s6m);
    sixsig->SetBinContent(ch+1,sig6);
    sixsig->SetBinError(ch+1,s6e);
    std::string plname="plots/auto-run68-baseRMS-gausFit-ch"+std::to_string(ch)+".png";
    c1.Print(plname.c_str());
   }
  c1.Clear();
  gStyle->SetOptStat(0);
  mean->Draw();
  c1.Print("plots/auto-run68-baseRMS-meanFit-vs-ch.png");

  c1.Clear();
  sig->Draw();
  c1.Print("plots/auto-run68-baseRMS-sigmaFit-vs-ch.png");

  c1.Clear();
  sixsig->Draw();
  c1.Print("plots/auto-run68-baseRMS-6sigmaCut-vs-ch.png");
  c1.Print("plots/auto-run68-baseRMS-6sigmaCut-vs-ch.C");
  return 0;
 }
