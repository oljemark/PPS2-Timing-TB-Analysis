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
  float htotal2[11]={1711, 2836, 7600, 2341, 987, 507, 311, 345, 287, 174, 79};
  float htotal3[11]={177, 727, 2222, 1195, 499, 259, 180, 212, 168, 79, 62};
  float htotal4[11]={12, 554, 2539, 1233, 466, 221, 139, 143, 132, 61, 37};
  float htotal5[11]={2752, 5178, 14305, 4186, 1574, 808, 546, 694, 552, 273, 183};
  float htotal6[11]={1448, 9221, 16619, 5015, 1949, 1121, 735, 883, 615, 322, 181};
  float htotal7[11]={96, 321, 1329, 850, 591, 311, 194, 244, 175, 107, 79};
  float htotal8[11]={14158, 1322, 1312, 841, 393, 223, 163, 161, 139, 81, 36};
  float htotal9[11]={43, 394, 1014, 667, 337, 204, 126, 172, 125, 45, 32};
  float hpmt0[11]={0, 70, 965, 857, 347, 180, 120, 162, 114, 70, 43};
  float hpmt1[11]={0, 415, 3537, 3263, 1404, 839, 548, 664, 493, 276, 125};
  float hpmt2[11]={0, 2631, 7327, 2110, 812, 446, 295, 336, 280, 174, 77};
  float hpmt3[11]={0, 633, 2049, 1069, 456, 249, 179, 207, 167, 79, 62};
  float hpmt4[11]={0, 468, 2355, 1039, 382, 190, 133, 142, 128, 61, 37};
  float hpmt5[11]={0, 4835, 13923, 4046, 1526, 795, 535, 686, 548, 273, 180};
  float hpmt6[11]={0, 8936, 16083, 4728, 1861, 1101, 719, 877, 610, 315, 179};
  float hpmt7[11]={0, 170, 1271, 776, 549, 299, 188, 239, 173, 107, 76};
  float hpmt8[11]={1, 670, 1151, 730, 368, 214, 161, 160, 132, 81, 35};
  float hpmt9[11]={0, 221, 893, 606, 314, 199, 122, 171, 125, 43, 32};


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
  c1.Print("plots/Run44-PMT-coincidence-efficiency-vs-amplitude-multipad-effH.png");
  c1.Clear();
  pure0->SetLineColor(kRed);
  pure0->Draw();
  pure3->SetLineColor(kGreen);
  pure3->Draw("same");
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
  c1.Print("plots/Run44-PMT-coincidence-efficiency-vs-amplitude-single-pad-effH.png");
 }
