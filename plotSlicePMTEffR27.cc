#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include <cmath>

void plotSlicePMTEffR27()
 {
  float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};
  float htotal0[11]={  61,   30,   51,   47,  15,   1,   0,   0,   0,   0,   0};
  float htotal1[11]={1151, 1380, 4774, 1955, 701, 381, 261, 333, 264, 145, 115};
  float htotal2[11]={ 325,  932, 3328, 1059, 351, 202, 125, 165, 120,  75,  40};
  float htotal3[11]={  93,  363, 1090,  480, 194,  96,  77,  84,  80,  48,  31};
  float htotal4[11]={  18,  303, 1156,  450, 155,  76,  62,  74,  48,  29,  14};
  float htotal5[11]={   2,  262, 2966, 1791, 702, 442, 272, 310, 245, 158,  98};
  float htotal6[11]={  18, 2981, 8375, 2649, 958, 519, 354, 458, 334, 192, 112};
  float htotal7[11]={  50,  129, 1019,  478, 249, 169, 106, 134,  95,  55,  35};
  float htotal8[11]={   2,   58,  767,  602, 355, 184, 114, 152, 101,  84,  54};
  float htotal9[11]={6587, 2604, 3864,  781, 328, 149,  86, 123, 106,  57,  37};
  float hpmt0[11]=  {   0,    0,    0,    0,   0,   0,   0,   0,   0,   0,   0};
  float hpmt1[11]=  {   0, 1324, 4563, 1787, 662, 372, 260, 330, 257, 143, 112};
  float hpmt2[11]=  {   0,  850, 3221,  970, 318, 192, 124, 164, 120,  74,  40};
  float hpmt3[11]=  {   0,  297, 1023,  451, 185,  94,  76,  82,  77,  47,  30};
  float hpmt4[11]=  {   0,  243, 1052,  377, 130,  74,  61,  71,  48,  28,  14};
  float hpmt5[11]=  {   0,  222, 2823, 1727, 683, 437, 271, 306, 241, 157,  97};
  float hpmt6[11]=  {   0, 2885, 8137, 2491, 933, 511, 349, 454, 330, 191, 110};
  float hpmt7[11]=  {   0,   74,  978,  440, 233, 144,  86, 121,  95,  55,  35};
  float hpmt8[11]=  {   0,   42,  708,  520, 310, 165,  93, 127, 101,  83,  52};
  float hpmt9[11]=  {   0, 2297, 3717,  690, 259, 119,  82, 122, 106,  57,  37};


  auto TW0=new TH1F("totWaves_0","Total waveforms (ch.0);Amplitude (mV)", 11, binEdges);
  auto TW1=new TH1F("totWaves_1","Total waveforms (ch.1);Amplitude (mV)", 11, binEdges);
  auto TW2=new TH1F("totWaves_2","Total waveforms (ch.2);Amplitude (mV)", 11, binEdges);
  auto TW3=new TH1F("totWaves_3","Total waveforms (ch.3);Amplitude (mV)", 11, binEdges);
  auto TW4=new TH1F("totWaves_4","Total waveforms (ch.4);Amplitude (mV)", 11, binEdges);
  auto TW5=new TH1F("totWaves_5","Total waveforms (ch.5);Amplitude (mV)", 11, binEdges);
  auto TW6=new TH1F("totWaves_6","Total waveforms (ch.6);Amplitude (mV)", 11, binEdges);
  auto TW7=new TH1F("totWaves_7","Total waveforms (ch.7);Amplitude (mV)", 11, binEdges);
  auto TW8=new TH1F("totWaves_8","Total waveforms (ch.8);Amplitude (mV)", 11, binEdges);
  auto TW9=new TH1F("totWaves_9","Total waveforms (ch.9);Amplitude (mV)", 11, binEdges);

  auto TW0coi=new TH1F("pmtCoincidentWaves_0","Waveforms with PMT coincidence (ch.0);Amplitude (mV)", 11, binEdges);
  auto TW1coi=new TH1F("pmtCoincidentWaves_1","Waveforms with PMT coincidence (ch.1);Amplitude (mV)", 11, binEdges);
  auto TW2coi=new TH1F("pmtCoincidentWaves_2","Waveforms with PMT coincidence (ch.2);Amplitude (mV)", 11, binEdges);
  auto TW3coi=new TH1F("pmtCoincidentWaves_3","Waveforms with PMT coincidence (ch.3);Amplitude (mV)", 11, binEdges);
  auto TW4coi=new TH1F("pmtCoincidentWaves_4","Waveforms with PMT coincidence (ch.4);Amplitude (mV)", 11, binEdges);
  auto TW5coi=new TH1F("pmtCoincidentWaves_5","Waveforms with PMT coincidence (ch.5);Amplitude (mV)", 11, binEdges);
  auto TW6coi=new TH1F("pmtCoincidentWaves_6","Waveforms with PMT coincidence (ch.6);Amplitude (mV)", 11, binEdges);
  auto TW7coi=new TH1F("pmtCoincidentWaves_7","Waveforms with PMT coincidence (ch.7);Amplitude (mV)", 11, binEdges);
  auto TW8coi=new TH1F("pmtCoincidentWaves_8","Waveforms with PMT coincidence (ch.8);Amplitude (mV)", 11, binEdges);
  auto TW9coi=new TH1F("pmtCoincidentWaves_9","Waveforms with PMT coincidence (ch.9);Amplitude (mV)", 11, binEdges);

  for (int i=0;i < 11;i++)
   {
    TW1->SetBinContent(i+1,htotal1[i]);
    TW6->SetBinContent(i+1,htotal6[i]);
    TW1coi->SetBinContent(i+1,hpmt1[i]);
    TW6coi->SetBinContent(i+1,hpmt6[i]);
    TW8->SetBinContent(i+1,htotal8[i]);
    TW8coi->SetBinContent(i+1,hpmt8[i]);
    TW2->SetBinContent(i+1,htotal2[i]);
    TW2coi->SetBinContent(i+1,hpmt2[i]);
    TW3->SetBinContent(i+1,htotal3[i]);
    TW3coi->SetBinContent(i+1,hpmt3[i]);
    TW4->SetBinContent(i+1,htotal4[i]);
    TW4coi->SetBinContent(i+1,hpmt4[i]);
    TW5->SetBinContent(i+1,htotal5[i]);
    TW5coi->SetBinContent(i+1,hpmt5[i]);
    TW7->SetBinContent(i+1,htotal7[i]);
    TW7coi->SetBinContent(i+1,hpmt7[i]);
    TW9->SetBinContent(i+1,htotal9[i]);
    TW9coi->SetBinContent(i+1,hpmt9[i]);
    TW0coi->SetBinContent(i+1,hpmt0[i]);
    TW0->SetBinContent(i+1,htotal0[i]);
   }
  TW1->Sumw2();
  TW6->Sumw2();
  TW1coi->Sumw2();
  TW6coi->Sumw2();
  TW8->Sumw2();
  TW8coi->Sumw2();
  TW0->Sumw2();
  TW0coi->Sumw2();
  TW2->Sumw2();
  TW2coi->Sumw2();
  TW3->Sumw2();
  TW3coi->Sumw2();
  TW4->Sumw2();
  TW4coi->Sumw2();
  TW5->Sumw2();
  TW5coi->Sumw2();
  TW7->Sumw2();
  TW7coi->Sumw2();
  TW9->Sumw2();
  TW9coi->Sumw2();
  auto pure1=new TEfficiency(*TW1coi,*TW1);
  auto pure6=new TEfficiency(*TW6coi,*TW6);
  auto pure8=new TEfficiency(*TW8coi,*TW8);
  auto pure0=new TEfficiency(*TW0coi,*TW0);
  auto pure2=new TEfficiency(*TW2coi,*TW2);
  auto pure3=new TEfficiency(*TW3coi,*TW3);
  auto pure4=new TEfficiency(*TW4coi,*TW4);
  auto pure5=new TEfficiency(*TW5coi,*TW5);
  auto pure7=new TEfficiency(*TW7coi,*TW7);
  auto pure9=new TEfficiency(*TW9coi,*TW9);
  TW6->SetLineColor(kRed);
  TW1->SetLineColor(kRed);
  TW8->SetLineColor(kRed);
  TW0->SetLineColor(kRed);
  TW2->SetLineColor(kRed);
  TW3->SetLineColor(kRed);
  TW4->SetLineColor(kRed);
  TW5->SetLineColor(kRed);
  TW7->SetLineColor(kRed);
  TW9->SetLineColor(kRed);

  TCanvas c1;
  gStyle->SetOptStat(0);
  TW1->Draw();
  TW1coi->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->AddEntry(TW1,"all amplitudes (ch.1)","l");
  leg->AddEntry(TW1coi,"PMT coincidence (ch.1)","l");
  leg->Draw();
  c1.Print("plots/Run27-amplitude-slices-PMT-eff-ch1.png");
  c1.Clear();
  TW1coi->Divide(TW1);

  TW9->Draw();
  TW9coi->Draw("same");
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(TW9,"all amplitudes (ch.9)","l");
  leg2->AddEntry(TW9coi,"PMT coincidence (ch.9)","l");
  leg2->Draw();
  c1.Print("plots/Run27-amplitude-slices-PMT-eff-ch9.png");
  c1.Clear();
  TW9coi->Divide(TW9);
  TW6->Draw();
  TW6coi->Draw("same");
  auto leg4=new TLegend(0.7,0.7,0.9,0.9);
  leg4->AddEntry(TW6,"all amplitudes (ch.6)","l");
  leg4->AddEntry(TW6coi,"PMT coincidence (ch.6)","l");
  leg4->Draw();
  c1.Print("plots/Run27-amplitude-slices-PMT-eff-ch6.png");
  c1.Clear();
//  TW6coi->Divide(TW6);
//  TW1coi->Draw();
  pure1->SetLineColor(kRed);
  pure1->Draw();
  pure2->SetLineColor(kGreen);
  pure2->Draw("same");
  pure5->SetLineColor(kCyan);
  pure5->Draw("same");
  pure6->SetLineColor(kBlack);
  pure6->Draw("same");
  auto leg3=new TLegend(0.7,0.3,0.9,0.6);
  leg3->SetHeader("Multipad ch.","C");
  leg3->AddEntry(pure1,"channel 1 purity","l");
  leg3->AddEntry(pure2,"channel 2 purity","l");
  leg3->AddEntry(pure5,"channel 5 purity","l");
  leg3->AddEntry(pure6,"channel 6 purity","l");
  leg3->Draw();
  c1.Print("plots/Run27-PMT-coincidence-efficiency-vs-amplitude-multipad-effH.png");
  c1.Clear();
  pure3->SetLineColor(kGreen);
  pure3->Draw();
  pure0->SetLineColor(kRed);
  pure0->Draw("same");
  pure4->SetLineColor(kCyan);
  pure4->Draw("same");
  pure7->SetLineColor(kBlack);
  pure7->Draw("same");
  pure8->SetLineColor(kYellow);
  pure8->Draw("same");
  pure9->SetLineColor(kBrown);
  pure9->Draw("same");
  auto leg5=new TLegend(0.7,0.3,0.9,0.6);
  leg5->SetHeader("Single pad ch.","C");
  leg5->AddEntry(pure0,"channel 0 purity","l");
  leg5->AddEntry(pure3,"channel 3 purity","l");
  leg5->AddEntry(pure4,"channel 4 purity","l");
  leg5->AddEntry(pure7,"channel 7 purity","l");
  leg5->AddEntry(pure8,"channel 8 purity","l");
  leg5->AddEntry(pure9,"channel 9 purity","l");
  leg5->Draw();
  c1.Print("plots/Run27-PMT-coincidence-efficiency-vs-amplitude-single-pad-effH.png");
 }
