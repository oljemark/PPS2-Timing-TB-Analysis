#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>

void plotNoisyBaseline()
 {
  const float waveForms[30]={36621, 13358, 113302, 35872, 12676, 47328,  22347, 33332,  20170, 11529,
                             23405, 13485,  11287, 22110,  5130, 51390, 173874, 74232, 111411, 63543,
                             16334, 17951,  34416, 35751, 23910, 26858,  11988, 14430,  11251, 26713};
  const float noisy[30]=    { 1302,   877,   2386,  1496,   974,  1331,   1235,  1300,  1039,    770,
                              1074,  1142,    778,  1333,   863,  1761,   1761,  1411,  1720,  1116,
                              1098,   834,   1212,  1002,   928,  3201,   2233,  1242,  2487,  1863};


  auto wave=new TH1F("totWaves","Total waveforms;Channel num.", 30, -0.5,29.5);
  auto noise=new TH1F("noiseBase","Noisy baseline waveforms (RMS > 6sig, ~2 mV);Channel num.", 30, -0.5,29.5);

  for (int i=0;i < 30;i++)
   {
    wave->SetBinContent(i+1,waveForms[i]);
    noise->SetBinContent(i+1,noisy[i]);
   }
  wave->Sumw2();
  noise->Sumw2();
  noise->SetLineColor(kRed);

  TCanvas c1;
  gStyle->SetOptStat(0);
  wave->GetYaxis()->SetRangeUser(0.,180000.);
  wave->Draw();
  noise->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->AddEntry(wave,"all waveforms","l");
  leg->AddEntry(noise,"baseline RMS above 6 sigma cut","l");
  leg->Draw();
  c1.Print("plots/Run68-baseline-noise-per-channel.png");
  c1.Clear();
  noise->Divide(wave);

  noise->Draw();
  auto leg2=new TLegend(0.1,0.8,0.3,0.9);
  leg2->AddEntry(noise,"noisy waveform fraction","l");
  leg2->Draw();
  c1.Print("plots/Run68-channel-baseline-noise-rates.png");
 }
