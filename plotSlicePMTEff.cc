#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include <cmath>

void plotSlicePMTEff()
 {
  float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};
  float htotal0[11]={3, 71, 1051, 1247, 639, 344, 190, 182, 114, 73, 43};
  float htotal1[11]={30, 464, 3793, 3638, 1515, 869, 559, 673, 500, 279, 128};
  float htotal6[11]={1448, 9221, 16619, 5015, 1949, 1121, 735, 883, 615, 322, 181};
  float htotal8[11]={14158, 1322, 1312, 841, 393, 223, 163, 161, 139, 81, 36};
  float hpmt0[11]={0, 70, 965, 857, 347, 180, 120, 162, 114, 70, 43};
  float hpmt1[11]={0, 415, 3537, 3263, 1404, 839, 548, 664, 493, 276, 125};
  float hpmt6[11]={0, 8936, 16083, 4728, 1861, 1101, 719, 877, 610, 315, 179};
  float hpmt8[11]={1, 670, 1151, 730, 368, 214, 161, 160, 132, 81, 35};


  auto TW0=new TH1F("totWaves_0","Total waveforms (ch.0);Amplitude (mV)", 11, binEdges);
  auto TW1=new TH1F("totWaves_1","Total waveforms (ch.1);Amplitude (mV)", 11, binEdges);
  auto TW6=new TH1F("totWaves_6","Total waveforms (ch.6);Amplitude (mV)", 11, binEdges);
  auto TW8=new TH1F("totWaves_8","Total waveforms (ch.8);Amplitude (mV)", 11, binEdges);

  auto TW0coi=new TH1F("pmtCoincidentWaves_0","Waveforms with PMT coincidence (ch.0);Amplitude (mV)", 11, binEdges);
  auto TW1coi=new TH1F("pmtCoincidentWaves_1","Waveforms with PMT coincidence (ch.1);Amplitude (mV)", 11, binEdges);
  auto TW6coi=new TH1F("pmtCoincidentWaves_6","Waveforms with PMT coincidence (ch.6);Amplitude (mV)", 11, binEdges);
  auto TW8coi=new TH1F("pmtCoincidentWaves_8","Waveforms with PMT coincidence (ch.8);Amplitude (mV)", 11, binEdges);

  for (int i=0;i < 11;i++)
   {
    TW1->SetBinContent(i+1,htotal1[i]);
    TW6->SetBinContent(i+1,htotal6[i]);
    TW1coi->SetBinContent(i+1,hpmt1[i]);
    TW6coi->SetBinContent(i+1,hpmt6[i]);
    TW8->SetBinContent(i+1,htotal8[i]);
    TW8coi->SetBinContent(i+1,hpmt8[i]);
    TW0->SetBinContent(i+1,htotal0[i]);
    TW0coi->SetBinContent(i+1,hpmt0[i]);
   }
  TW1->Sumw2();
  TW6->Sumw2();
  TW1coi->Sumw2();
  TW6coi->Sumw2();
  TW8->Sumw2();
  TW8coi->Sumw2();
  TW0->Sumw2();
  TW0coi->Sumw2();
  auto pure1=new TEfficiency(*TW1coi,*TW1);
  auto pure6=new TEfficiency(*TW6coi,*TW6);
  auto pure8=new TEfficiency(*TW8coi,*TW8);
  auto pure0=new TEfficiency(*TW0coi,*TW0);
  TW6->SetLineColor(kRed);
  TW1->SetLineColor(kRed);
  TW8->SetLineColor(kRed);
  TW0->SetLineColor(kRed);

  TCanvas c1;
  gStyle->SetOptStat(0);
  TW1->Draw();
  TW1coi->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->AddEntry(TW1,"all amplitudes (ch.1)","l");
  leg->AddEntry(TW1coi,"PMT coincidence (ch.1)","l");
  leg->Draw();
  c1.Print("plots/Run44-amplitude-slices-PMT-eff-ch1.png");
  c1.Clear();
  TW1coi->Divide(TW1);

  TW8->Draw();
  TW8coi->Draw("same");
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(TW8,"all amplitudes (ch.8)","l");
  leg2->AddEntry(TW8coi,"PMT coincidence (ch.8)","l");
  leg2->Draw();
  c1.Print("plots/Run44-amplitude-slices-PMT-eff-ch8.png");
  c1.Clear();
  TW8coi->Divide(TW8);
  TW6->Draw();
  TW6coi->Draw("same");
  auto leg4=new TLegend(0.7,0.7,0.9,0.9);
  leg4->AddEntry(TW6,"all amplitudes (ch.6)","l");
  leg4->AddEntry(TW6coi,"PMT coincidence (ch.6)","l");
  leg4->Draw();
  c1.Print("plots/Run44-amplitude-slices-PMT-eff-ch6.png");
  c1.Clear();
//  TW6coi->Divide(TW6);
//  TW1coi->Draw();
  pure8->SetLineColor(kRed);
  pure8->Draw();
  pure6->SetLineColor(kGreen);
  pure6->Draw("same");
  pure1->SetLineColor(kCyan);
  pure1->Draw("same");
  pure0->SetLineColor(kBlack);
  pure0->Draw("same");
  auto leg3=new TLegend(0.7,0.3,0.9,0.6);
  leg3->AddEntry(pure0,"channel 0 purity","l");
  leg3->AddEntry(pure1,"channel 1 purity","l");
  leg3->AddEntry(pure6,"channel 6 purity","l");
  leg3->AddEntry(pure8,"channel 8 purity","l");
  leg3->Draw();
  c1.Print("plots/Run44-PMT-coincidence-efficiency-vs-amplitude-ch0-ch1-ch6-ch8-effH.png");
 }
