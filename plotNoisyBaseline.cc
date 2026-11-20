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
  const float waves44[10]=  { 3960, 12454,  17243,  5787,  5541, 31112,  38234,  4345, 19380,  3160};
  const float noisy44[10]=  { 1013,   827,    787,   521,   528,   979,   1379,   427,   784,   349};

  auto wave=new TH1F("totWaves","Total waveforms, run 68;Channel num.", 30, -0.5,29.5);
  auto noise=new TH1F("noiseBase","Noisy baseline waveforms (RMS > 6sig), run 68;Channel num.", 30, -0.5,29.5);
  auto wave44=new TH1F("totWaves44","Total waveforms (Run 44, LD_04);Channel num.", 10, -0.5,9.5);
  auto noise44=new TH1F("noiseBase44","Noisy baseline waveforms (RMS > 6sig), run 44 LD_04;Channel num.", 10, -0.5,9.5);
  auto wave68_4=new TH1F("totWaves68_4","Total waveforms (Run 68, LD_04);Channel order as in run 44", 10, -0.5,9.5);
  auto noise68_4=new TH1F("noiseBase68_4","Noisy baseline waveforms (RMS > 6sig), run 68 LD_04;Channel order as in run 44", 10, -0.5,9.5);
  auto wave68_7=new TH1F("totWaves68_7","Total waveforms (Run 68, LD_07);Channel num. (in run 44)", 10, -0.5,9.5);
  auto noise68_7=new TH1F("noiseBase68_7","Noisy baseline waveforms (RMS > 6sig), run 68 LD_07;Channel num. (in run 44)", 10, -0.5,9.5);
  auto wave68_8=new TH1F("totWaves68_8","Total waveforms (Run 68, LD_08);Channel num. (in run 44)", 10, -0.5,9.5);
  auto noise68_8=new TH1F("noiseBase68_8","Noisy baseline waveforms (RMS > 6sig), run 68 LD_08;Channel num. (in run 44)", 10, -0.5,9.5);

  for (int i=0;i < 30;i++)
   {
    wave->SetBinContent(i+1,waveForms[i]);
    noise->SetBinContent(i+1,noisy[i]);
   }
  {
    wave68_4->SetBinContent(0+1,waveForms[4]);
    noise68_4->SetBinContent(0+1,noisy[4]);
    wave68_4->SetBinContent(1+1,waveForms[5]);
    noise68_4->SetBinContent(1+1,noisy[5]);
    wave68_7->SetBinContent(0+1,waveForms[6]);
    noise68_7->SetBinContent(0+1,noisy[6]);
    wave68_7->SetBinContent(1+1,waveForms[7]);
    noise68_7->SetBinContent(1+1,noisy[7]);
    wave68_8->SetBinContent(0+1,waveForms[1]);
    noise68_8->SetBinContent(0+1,noisy[1]);
    wave68_8->SetBinContent(1+1,waveForms[2]);
    noise68_8->SetBinContent(1+1,noisy[2]);

    wave68_4->SetBinContent(2+1,waveForms[0]);
    noise68_4->SetBinContent(2+1,noisy[0]);
    wave68_4->SetBinContent(3+1,waveForms[10]);
    noise68_4->SetBinContent(3+1,noisy[10]);
    wave68_7->SetBinContent(2+1,waveForms[12]);
    noise68_7->SetBinContent(2+1,noisy[12]);
    wave68_7->SetBinContent(3+1,waveForms[13]);
    noise68_7->SetBinContent(3+1,noisy[13]);
    wave68_8->SetBinContent(2+1,waveForms[3]);
    noise68_8->SetBinContent(2+1,noisy[3]);
    wave68_8->SetBinContent(3+1,waveForms[8]);
    noise68_8->SetBinContent(3+1,noisy[8]);

    wave68_4->SetBinContent(4+1,waveForms[11]);
    noise68_4->SetBinContent(4+1,noisy[11]);
    wave68_4->SetBinContent(5+1,waveForms[18]);
    noise68_4->SetBinContent(5+1,noisy[18]);
    wave68_7->SetBinContent(4+1,waveForms[14]);
    noise68_7->SetBinContent(4+1,noisy[14]);
    wave68_7->SetBinContent(5+1,waveForms[21]);
    noise68_7->SetBinContent(5+1,noisy[21]);
    wave68_8->SetBinContent(4+1,waveForms[9]);
    noise68_8->SetBinContent(4+1,noisy[9]);
    wave68_8->SetBinContent(5+1,waveForms[15]);
    noise68_8->SetBinContent(5+1,noisy[15]);

    wave68_4->SetBinContent(6+1,waveForms[19]);
    noise68_4->SetBinContent(6+1,noisy[19]);
    wave68_4->SetBinContent(7+1,waveForms[20]);
    noise68_4->SetBinContent(7+1,noisy[20]);
    wave68_7->SetBinContent(6+1,waveForms[22]);
    noise68_7->SetBinContent(6+1,noisy[22]);
    wave68_7->SetBinContent(7+1,waveForms[27]);
    noise68_7->SetBinContent(7+1,noisy[27]);
    wave68_8->SetBinContent(6+1,waveForms[16]);
    noise68_8->SetBinContent(6+1,noisy[16]);
    wave68_8->SetBinContent(7+1,waveForms[17]);
    noise68_8->SetBinContent(7+1,noisy[17]);

    wave68_4->SetBinContent(8+1,waveForms[25]);
    noise68_4->SetBinContent(8+1,noisy[25]);
    wave68_4->SetBinContent(9+1,waveForms[26]);
    noise68_4->SetBinContent(9+1,noisy[26]);
    wave68_7->SetBinContent(8+1,waveForms[28]);
    noise68_7->SetBinContent(8+1,noisy[28]);
    wave68_7->SetBinContent(9+1,waveForms[29]);
    noise68_7->SetBinContent(9+1,noisy[29]);
    wave68_8->SetBinContent(8+1,waveForms[23]);
    noise68_8->SetBinContent(8+1,noisy[23]);
    wave68_8->SetBinContent(9+1,waveForms[24]);
    noise68_8->SetBinContent(9+1,noisy[24]);

  }
  wave68_4->Sumw2();
  noise68_4->Sumw2();
  noise68_4->SetLineColor(kBlue);

  wave68_7->Sumw2();
  noise68_7->Sumw2();
  noise68_7->SetLineColor(kMagenta);

  wave68_8->Sumw2();
  noise68_8->Sumw2();
  noise68_8->SetLineColor(kCyan);

  wave->Sumw2();
  noise->Sumw2();
  noise->SetLineColor(kRed);

  for (int i=0;i < 10;i++)
   {
    wave44->SetBinContent(i+1,waves44[i]);
    noise44->SetBinContent(i+1,noisy44[i]);
   }
  wave44->Sumw2();
  noise44->Sumw2();
  noise44->SetLineColor(kRed);

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
  noise44->Divide(wave44);
  noise68_4->Divide(wave68_4);
  noise68_7->Divide(wave68_7);
  noise68_8->Divide(wave68_8);

  noise->Draw();
  auto leg2=new TLegend(0.1,0.8,0.3,0.9);
  leg2->AddEntry(noise,"noisy waveform fraction","l");
  leg2->Draw();
  c1.Print("plots/Run68-channel-baseline-noise-rates.png");
  c1.Print("plots/Run68-channel-baseline-noise-rates.C");
  c1.Clear();

  noise68_4->GetYaxis()->SetRangeUser(0.0,0.27);
  noise68_4->Draw();
  noise44->Draw("same");
  noise68_7->Draw("same");
  noise68_8->Draw("same");
  auto leg4=new TLegend(0.1,0.5,0.3,0.7);
  leg4->AddEntry(noise44,"LD_04 (Run 44)","l");
  leg4->AddEntry(noise68_4,"LD_04 (Run 68)","l");
  leg4->AddEntry(noise68_7,"LD_07 (Run 68)","l");
  leg4->AddEntry(noise68_8,"LD_08 (Run 68)","l");
  leg4->Draw();
  c1.Print("plots/Run68-vs-Run44-channel-baseline-noise-rates.png");
 }
