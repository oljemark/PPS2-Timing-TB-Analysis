#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include <cmath>

int preTWresolution()
 {
  auto rd=TFile::Open("triple-run068-lgad04-07-08-3GeV-40pct-v1.root","READ");
  if (!rd->IsOpen())
   return 1;
  auto mean=new TH1F("PMTtimeOffsetMean","Mean PMT time delay vs LGAD channel;Channel (number);Fitted mean (ps)",30,-0.5,29.5);
  auto sig=new TH1F("PMTtimeOffsetResolution","PMT time delay vs LGAD channel, fitted resolution;Channel (number);Fitted sigma (ps)",30,-0.5,29.5);
  TCanvas c1;
  gStyle->SetOptFit(1);
  for (uint8_t ch=0; ch < 30; ch++)
   {
    std::string hname="timeDiffPMT_ch32_"+std::to_string(ch)+";1";
    c1.Clear();
    auto rms=rd->Get<TH1F>(hname.c_str());
    auto resp=rms->Fit("gaus","S");
    mean->SetBinContent(ch+1,resp->Parameter(1));
    mean->SetBinError(ch+1,resp->ParError(1));
    sig->SetBinContent(ch+1,resp->Parameter(2));
    sig->SetBinError(ch+1,resp->ParError(2));
    std::string plname="plots/auto-run68-PMT-delay-gausFit-ch"+std::to_string(ch)+".png";
    c1.Print(plname.c_str());
   }
  c1.Clear();
  gStyle->SetOptStat(0);
  mean->Draw();
  c1.Print("plots/auto-run68-pmtDelay-meanFit-vs-ch.png");

  c1.Clear();
  sig->Draw();
  c1.Print("plots/auto-run68-pmtDelay-resolution-sigmaFit-vs-ch.png");

  return 0;
 }
