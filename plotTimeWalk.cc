#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>

void plotTimeWalk()
 {
  float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};
  float mean_ch1_40[11]={0, 682.24, 684.268, 695.72, 710.555, 720.937, 724.642, 729.911, 738.118, 736.52, 742.36};
  float meanE_ch1_40[11]={0, 4.99169, 1.61315, 1.91561, 3.51138, 4.96849, 6.62041, 5.13494, 8.30385, 7.46909, 12.5606};
  float rms_ch1_40[11]={0, 87.282, 95.1748, 107.809, 123.211, 131.671, 130.956, 117.77, 140.212, 112.573, 114.249};
  float rmsE_ch1_40[11]={0, 3.67151 , 1.10091, 1.28148, 2.25149, 3.31707, 5.1126, 3.34739, 7.51714, 6.16275, 11.0715};

  float mean_ch6_40[11]={0, 525.514, 541.634, 558.323, 571.407, 573.857, 586.542, 591.734, 590.836, 592.493, 597.086};
  float meanE_ch6_40[11]={0, 0.951367, 0.750533, 1.6382, 2.95744, 4.18269, 5.46047, 5.03878, 5.10588, 7.37598, 9.53549};
  float rms_ch6_40[11]={0, 88.8072, 94.8971, 110.998, 121.797, 127.323, 131.647, 132.771, 114.934, 115.699, 110.008};
  float rmsE_ch6_40[11]={0, 0.722861, 0.527734, 1.04295, 1.77913, 2.67377, 3.65539, 3.24536, 3.34186, 5.59873, 8.38311};

  auto TW1=new TH1F("timeWalkCorr_1","Time walk correction? Offset from highest amplitude (ch.1);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW6=new TH1F("timeWalkCorr_6","Time walk correction? Offset from highest amplitude (ch.6);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);

  auto mean40_1=new TH1F("mean40_1","Fitted mean time difference vs PMT (ch.1);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_1=new TH1F("rms40_1","Fitted RMS of time difference vs PMT (ch.1);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_6=new TH1F("mean40_6","Fitted mean time difference vs PMT (ch.6);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_6=new TH1F("rms40_6","Fitted RMS of time difference vs PMT (ch.6);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);

  for (int i=1;i < 11;i++)
   {
    float em1=meanE_ch1_40[i];
    float em1b=meanE_ch1_40[10];
    mean40_1->SetBinContent(i+1,mean_ch1_40[i]);
    mean40_1->SetBinError(i+1,meanE_ch1_40[i]);
    TW1->SetBinContent(i+1,mean_ch1_40[10]-mean_ch1_40[i]);
    TW1->SetBinError(i+1,sqrt(em1*em1+em1b*em1b));

    rms40_1->SetBinContent(i+1,rms_ch1_40[i]);
    rms40_1->SetBinError(i+1,rmsE_ch1_40[i]);
    mean40_6->SetBinContent(i+1,mean_ch6_40[i]);
    mean40_6->SetBinError(i+1,meanE_ch6_40[i]);
    rms40_6->SetBinContent(i+1,rms_ch6_40[i]);
    rms40_6->SetBinError(i+1,rmsE_ch6_40[i]);
    float em6=meanE_ch6_40[i];
    float em6b=meanE_ch6_40[10];

    TW6->SetBinContent(i+1,mean_ch6_40[10]-mean_ch6_40[i]);
    TW6->SetBinError(i+1,sqrt(em6*em6+em6b*em6b));

   }
  mean40_6->SetLineColor(kRed);
  rms40_6->SetLineColor(kRed);
  TW6->SetLineColor(kRed);

  TCanvas c1;
  gStyle->SetOptStat(0);
  mean40_1->GetYaxis()->SetRangeUser(500.,750.);
  rms40_1->GetYaxis()->SetRangeUser(80.,150.);
  TW1->GetYaxis()->SetRangeUser(-10.,100.);
  mean40_1->Draw();
  mean40_6->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->AddEntry(mean40_1,"channel 1","l");
  leg->AddEntry(mean40_6,"channel 6","l");
  leg->Draw();
  c1.Print("plots/Run44-timeWalk-meanTime-vs-amplitude-ch1-ch6.png");
  c1.Clear();
  rms40_1->Draw();
  rms40_6->Draw("same");
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(rms40_1,"channel 1","l");
  leg2->AddEntry(rms40_6,"channel 6","l");
  leg2->Draw();
  c1.Print("plots/Run44-timeWalk-rmsTime-vs-amplitude-ch1-ch6.png");
  c1.Clear();
  TW1->Draw();
  TW6->Draw("same");
  auto leg3=new TLegend(0.7,0.7,0.9,0.9);
  leg3->AddEntry(TW1,"channel 1","l");
  leg3->AddEntry(TW6,"channel 6","l");
  leg3->Draw();
  c1.Print("plots/Run44-timeWalk-correction-vs-amplitude-ch1-ch6.png");
  c1.Clear();
  gStyle->SetOptFit(1);
  auto walk=new TF1("walk","[0]/TMath::Sqrt(x+[1])+[2]",5.,150.);
//  walk->SetParLimits(0,0.05,1000.);
//  walk->SetParLimits(1,-10.,10.);
//  walk->SetParLimits(2,-40.,40.);
  TW1->Fit("walk");
  c1.Print("plots/Run44-timeWalk-correction-fit5-ch1.png");
  TW6->Fit("walk");
  c1.Print("plots/Run44-timeWalk-correction-fit5-ch6.png");
  mean40_6->Fit("walk");
  c1.Print("plots/Run44-timeWalk-fit-ch6.png");
  mean40_1->Fit("walk");
  c1.Print("plots/Run44-timeWalk-fit-ch1.png");

 }
