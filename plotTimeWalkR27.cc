#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>

void plotTimeWalkR27()
 {
  float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};
  float mean_ch1_40[11]= {0, 688.225, 693.279, 699.977,  715.09, 688.058, 701.465, 742.508, 746.156, 748.851, 736.826};// done
  float meanE_ch1_40[11]={0, 2.63908, 1.32547,  2.7388, 4.98339, 9.41366, 16.3458, 8.69573, 8.76945, 17.6661, 12.4682};
  float rms_ch1_40[11]=  {0,  88.478, 88.1267, 110.085, 120.375, 140.302, 182.762, 142.381, 124.263, 149.907, 106.61};
  float rmsE_ch1_40[11]= {0, 1.97122, 1.00069, 1.90897, 3.44649, 7.49073, 17.4952, 6.89681, 6.92076, 17.3593, 12.2665};

  float mean_ch6_40[11]= {0, 517.563, 542.556, 559.945, 565.201, 568.516, 578.773, 600.735, 600.715, 605.962, 586.811};
  float meanE_ch6_40[11]={0, 1.62593, 1.00836, 2.22599, 4.34803, 6.47838, 8.22098, 6.78684, 8.07263,  9.7785, 10.6829};
  float rms_ch6_40[11]=  {0, 85.8773, 90.4701, 108.141,  124.31, 127.507,  129.99, 127.152, 128.042, 119.477, 98.9806};
  float rmsE_ch6_40[11]= {0, 1.24274, 0.70544, 1.52191,  2.8284, 4.60643, 5.53903, 4.86794, 6.19673, 7.60753, 10.2321};

  float mean_ch2_40[11]= {0, 585.367, 601.698, 615.977, 644.324,  633.86, 652.041, 654.124, 674.583, 662.557, 652.358};
  float meanE_ch2_40[11]={0, 2.8971,  1.50918, 3.57282, 8.59045, 13.1414,  25.934, 33.2979, 16.6181, 20.7264, 82.0939};
  float rms_ch2_40[11]=  {0, 80.8154, 84.9049, 104.235, 131.892, 148.078, 182.496, 219.577, 139.804, 130.658, 209.159};
  float rmsE_ch2_40[11]= {0, 2.26485, 1.12164, 2.48066, 7.25494, 11.2277, 29.4752, 33.0834, 17.4173, 25.7329, 122.328};

  float mean_ch5_40[11]= {0, 492.902, 513.225,  520.73, 538.774, 543.176, 553.922, 557.735, 569.637, 576.219, 561.685};
  float meanE_ch5_40[11]={0, 6.18452, 1.67782, 2.62206, 5.24439, 7.43951, 13.7394,  10.453, 9.57642, 11.2766, 17.0577};
  float rms_ch5_40[11]=  {0,  84.289, 88.2046, 106.553, 125.453, 136.581,  159.54, 144.494, 126.987, 120.839, 123.335};
  float rmsE_ch5_40[11]= {0, 4.87968, 1.21527, 1.81194,  3.5791, 5.48281, 11.1958,  8.2727, 8.70975, 10.2285, 17.2665};

  auto TW1=new TH1F("timeWalkCorr_1","Time walk correction? Offset from highest amplitude (ch.1);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW2=new TH1F("timeWalkCorr_2","Time walk correction? Offset from highest amplitude (ch.2);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW5=new TH1F("timeWalkCorr_5","Time walk correction? Offset from highest amplitude (ch.5);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW6=new TH1F("timeWalkCorr_6","Time walk correction? Offset from highest amplitude (ch.6);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);

  auto mean40_1=new TH1F("mean40_1","Fitted mean time difference vs PMT (ch.1);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_1=new TH1F("rms40_1","Fitted RMS of time difference vs PMT (ch.1);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_2=new TH1F("mean40_2","Fitted mean time difference vs PMT (ch.2);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_2=new TH1F("rms40_2","Fitted RMS of time difference vs PMT (ch.2);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_5=new TH1F("mean40_5","Fitted mean time difference vs PMT (ch.5);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_5=new TH1F("rms40_5","Fitted RMS of time difference vs PMT (ch.5);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_6=new TH1F("mean40_6","Fitted mean time difference vs PMT (ch.6);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_6=new TH1F("rms40_6","Fitted RMS of time difference vs PMT (ch.6);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);

  for (int i=1;i < 11;i++)
   {
    float em1=meanE_ch1_40[i];
    float em1b=meanE_ch1_40[8];
    mean40_1->SetBinContent(i+1,mean_ch1_40[i]);
    mean40_1->SetBinError(i+1,meanE_ch1_40[i]);
    TW1->SetBinContent(i+1,mean_ch1_40[8]-mean_ch1_40[i]);
    TW1->SetBinError(i+1,sqrt(em1*em1+em1b*em1b));
    rms40_1->SetBinContent(i+1,rms_ch1_40[i]);
    rms40_1->SetBinError(i+1,rmsE_ch1_40[i]);

    float em2=meanE_ch2_40[i];
    float em2b=meanE_ch2_40[8];
    mean40_2->SetBinContent(i+1,mean_ch2_40[i]);
    mean40_2->SetBinError(i+1,meanE_ch2_40[i]);
    TW2->SetBinContent(i+1,mean_ch2_40[8]-mean_ch2_40[i]);
    TW2->SetBinError(i+1,sqrt(em2*em2+em2b*em2b));
    rms40_2->SetBinContent(i+1,rms_ch2_40[i]);
    rms40_2->SetBinError(i+1,rmsE_ch2_40[i]);

    float em5=meanE_ch5_40[i];
    float em5b=meanE_ch5_40[8];
    mean40_5->SetBinContent(i+1,mean_ch5_40[i]);
    mean40_5->SetBinError(i+1,meanE_ch5_40[i]);
    TW5->SetBinContent(i+1,mean_ch5_40[8]-mean_ch5_40[i]);
    TW5->SetBinError(i+1,sqrt(em5*em5+em5b*em5b));
    rms40_5->SetBinContent(i+1,rms_ch5_40[i]);
    rms40_5->SetBinError(i+1,rmsE_ch5_40[i]);

    mean40_6->SetBinContent(i+1,mean_ch6_40[i]);
    mean40_6->SetBinError(i+1,meanE_ch6_40[i]);
    rms40_6->SetBinContent(i+1,rms_ch6_40[i]);
    rms40_6->SetBinError(i+1,rmsE_ch6_40[i]);
    float em6=meanE_ch6_40[i];
    float em6b=meanE_ch6_40[8];

    TW6->SetBinContent(i+1,mean_ch6_40[8]-mean_ch6_40[i]);
    TW6->SetBinError(i+1,sqrt(em6*em6+em6b*em6b));

   }
  mean40_6->SetLineColor(kRed);
  rms40_6->SetLineColor(kRed);
  TW6->SetLineColor(kRed);

  mean40_5->SetLineColor(kGreen);
  rms40_5->SetLineColor(kGreen);
  TW5->SetLineColor(kGreen);

  mean40_2->SetLineColor(kBlue);
  rms40_2->SetLineColor(kBlue);
  TW2->SetLineColor(kBlue);

  mean40_1->SetLineColor(kCyan);
  rms40_1->SetLineColor(kCyan);
  TW1->SetLineColor(kCyan);

  TCanvas c1;
  gStyle->SetOptStat(0);
  mean40_1->GetYaxis()->SetRangeUser(450.,800.);
  rms40_1->GetYaxis()->SetRangeUser(70.,230.);
  TW1->GetYaxis()->SetRangeUser(-40.,170.);
  mean40_1->Draw();
  mean40_2->Draw("same");
  mean40_5->Draw("same");
  mean40_6->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->SetHeader("Multipad ch.","C");
  leg->AddEntry(mean40_1,"channel 1","l");
  leg->AddEntry(mean40_2,"channel 2","l");
  leg->AddEntry(mean40_5,"channel 5","l");
  leg->AddEntry(mean40_6,"channel 6","l");
  leg->Draw();
  c1.Print("plots/Run27-timeWalk-meanTime-vs-amplitude-multiPads-1-2-5-6.png");
  c1.Clear();
  rms40_1->Draw();
  rms40_2->Draw("same");
  rms40_5->Draw("same");
  rms40_6->Draw("same");
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->SetHeader("Multipad ch.","C");
  leg2->AddEntry(rms40_1,"channel 1","l");
  leg2->AddEntry(rms40_2,"channel 2","l");
  leg2->AddEntry(rms40_5,"channel 5","l");
  leg2->AddEntry(rms40_6,"channel 6","l");
  leg2->Draw();
  c1.Print("plots/Run27-timeWalk-rmsTime-vs-amplitude-multiPads-1-2-5-6.png");
  c1.Clear();
  TW1->Draw();
  TW2->Draw("same");
  TW5->Draw("same");
  TW6->Draw("same");
  auto leg3=new TLegend(0.7,0.7,0.9,0.9);
  leg3->SetHeader("Multipad ch.","C");
  leg3->AddEntry(TW1,"channel 1","l");
  leg3->AddEntry(TW2,"channel 2","l");
  leg3->AddEntry(TW5,"channel 5","l");
  leg3->AddEntry(TW6,"channel 6","l");
  leg3->Draw();
  c1.Print("plots/Run27-timeWalk-correction-vs-amplitude-multiPads-1-2-5-6-45mV-ref.png");
  c1.Clear();
  gStyle->SetOptFit(1);
  auto walk=new TF1("walk","[0]/TMath::Sqrt(x+[1])+[2]",5.,150.);
  auto fex=new TF1("fex","expo(0)",5.,150.);
//  walk->SetParLimits(0,0.05,1000.);
//  walk->SetParLimits(1,-10.,10.);
//  walk->SetParLimits(2,-40.,40.);
  TW1->GetYaxis()->SetRangeUser(-20.,120.);
  TW1->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch1-45mV-ref.png");
  TW2->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch2-45mV-ref.png");
  TW5->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch5-45mV-ref.png");
  TW6->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch6-45mV-ref.png");
  mean40_6->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-ch6.png");
  mean40_5->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-ch5.png");
  mean40_2->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-ch2.png");
  mean40_1->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-ch1.png");

 }
