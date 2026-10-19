#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>

void plotNoisyBaseline()
 {
  const float waveForms[10]={3960, 12454, 17243, 5787, 5541, 31112, 38234, 4345, 19380, 3160};
  const float noisy[10]={992, 656, 661, 371, 479, 522, 814, 297, 371, 265};


  auto wave=new TH1F("totWaves","Total waveforms;Channel num.", 10, -0.5,9.5);
  auto noise=new TH1F("noiseBase","Noisy baseline waveforms (RMS > 2mV);Channel num.", 10, -0.5,9.5);

  for (int i=0;i < 10;i++)
   {
    wave->SetBinContent(i+1,waveForms[i]);
    noise->SetBinContent(i+1,noisy[i]);
   }
  wave->Sumw2();
  noise->Sumw2();
  noise->SetLineColor(kRed);

  TCanvas c1;
  gStyle->SetOptStat(0);
  wave->GetYaxis()->SetRangeUser(0.,40000.);
  wave->Draw();
  noise->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->AddEntry(wave,"all waveforms","l");
  leg->AddEntry(noise,"baseline RMS above 2 mV","l");
  leg->Draw();
  c1.Print("plots/Run44-baseline-noise-per-channel.png");
  c1.Clear();
  noise->Divide(wave);

  noise->Draw();
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(noise,"noisy waveform fraction","l");
  leg2->Draw();
  c1.Print("plots/Run44-channel-baseline-noise-rates.png");
 }
