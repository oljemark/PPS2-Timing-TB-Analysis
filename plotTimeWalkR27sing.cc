#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <cmath>

void plotTimeWalkR27sing()
 {
  float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};
  float mean_ch3_40[11]= {0, 429.541, 452.363, 461.755, 500.445, 509.691, 495.634, 507.228, 530.709, 555.905, 507.584};
  float meanE_ch3_40[11]={0, 4.90928, 2.57641, 5.56214, 12.3607, 31.9084, 41.6388, 25.5822, 21.8554, 26.2596, 40.0164};
  float rms_ch3_40[11]=  {0, 77.7668, 81.2314, 110.082, 136.994, 183.425, 208.413, 176.299, 141.546, 143.137, 153.115};
  float rmsE_ch3_40[11]= {0, 3.94849, 1.95432, 4.25905, 11.4298, 35.3253, 55.9743, 21.4794,  26.836, 26.7362, 43.2612};

  float mean_ch8_40[11]= {0, 169.507, 156.745, 172.927, 176.041, 201.177, 239.468, 192.912, 223.834, 251.549, 217.746};
  float meanE_ch8_40[11]={0, 15.8337, 2.72828, 4.79869, 9.59047, 14.2528,  19.709, 15.3873, 15.8685, 17.5276, 23.8319};
  float rms_ch8_40[11]=  {0, 72.2584, 71.3172, 101.711, 127.982, 136.788, 146.511, 149.389, 128.446, 134.988, 139.985};
  float rmsE_ch8_40[11]= {0, 25.1966, 1.90396,  3.7315, 7.75823,  12.824, 17.1852, 12.5102, 12.5422, 19.6482, 22.6879};

  float mean_ch4_40[11]= {0, 540.391, 555.216, 577.576, 582.033, 601.575, 658.534, 644.236,  654.63, 607.236, 682.568};
  float meanE_ch4_40[11]={0, 5.48447, 2.75354, 6.66718,   17.81, 19.4261, 45.9743, 21.8541, 20.7407, 43.0279, 35.0409};
  float rms_ch4_40[11]=  {0, 77.4236, 85.2155, 114.837, 156.262, 145.594, 186.682, 153.275, 128.518, 155.348, 107.577};
  float rmsE_ch4_40[11]= {0, 5.04347, 2.09296, 5.32838, 18.8366, 17.5483, 40.4318, 23.0977, 15.8584, 53.4749, 29.9006};

  float mean_ch7_40[11]= {0,  333.33, 342.903, 362.517, 371.993, 377.914, 402.783, 400.575, 402.386, 389.505, 432.829};
  float meanE_ch7_40[11]={0, 14.3321, 2.40049, 5.60038, 10.0426, 13.3915, 18.4564, 19.4152, 18.2294, 41.1197, 20.8132};
  float rms_ch7_40[11]=  {0, 86.9847, 73.6531, 105.953, 130.243, 137.235, 139.816, 158.197, 134.231,  202.51, 119.151};
  float rmsE_ch7_40[11]= {0, 16.0915, 1.82367, 5.02619,  8.4563, 9.41008, 18.0782, 21.8471, 16.4496,  55.852, 16.0307};

  float mean_ch9_40[11]= {0, 136.299,   143.22, 161.191, 161.604, 197.516, 199.088, 183.012,  218.51, 210.814, 176.162};
  float meanE_ch9_40[11]={0,  1.5993,  1.29755, 4.61995, 11.0701, 25.8868, 22.2738, 22.3288, 17.7146, 17.9893, 36.5973};
  float rms_ch9_40[11]=  {0, 75.2555,  78.3463, 110.478, 147.259, 183.983, 156.546, 162.882, 147.706, 111.307, 151.902};
  float rmsE_ch9_40[11]= {0, 1.29127, 0.984751, 3.51927, 11.3807,  26.431, 18.7127, 20.3631, 15.0323, 17.8424, 33.7214};

  auto TW3=new TH1F("timeWalkCorr_3","Time walk correction? Offset from highest amplitude (ch.3);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW4=new TH1F("timeWalkCorr_4","Time walk correction? Offset from highest amplitude (ch.4);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW7=new TH1F("timeWalkCorr_7","Time walk correction? Offset from highest amplitude (ch.7);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW8=new TH1F("timeWalkCorr_8","Time walk correction? Offset from highest amplitude (ch.8);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);
  auto TW9=new TH1F("timeWalkCorr_9","Time walk correction? Offset from highest amplitude (ch.9);Amplitude (mV);Time walk correction (difference vs largest bin, ps)", 11, binEdges);

  auto mean40_3=new TH1F("mean40_3","Fitted mean time difference vs PMT (ch.3);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_3=new TH1F("rms40_3","Fitted RMS of time difference vs PMT (ch.3);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_4=new TH1F("mean40_4","Fitted mean time difference vs PMT (ch.4);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_4=new TH1F("rms40_4","Fitted RMS of time difference vs PMT (ch.4);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_7=new TH1F("mean40_7","Fitted mean time difference vs PMT (ch.7);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_7=new TH1F("rms40_7","Fitted RMS of time difference vs PMT (ch.7);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_8=new TH1F("mean40_8","Fitted mean time difference vs PMT (ch.8);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_8=new TH1F("rms40_8","Fitted RMS of time difference vs PMT (ch.8);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);
  auto mean40_9=new TH1F("mean40_9","Fitted mean time difference vs PMT (ch.9);Amplitude (mV);Time difference LGAD - PMT (ps)", 11, binEdges);
  auto rms40_9=new TH1F("rms40_9","Fitted RMS of time difference vs PMT (ch.9);Amplitude (mV);RMS time difference LGAD - PMT (ps)", 11, binEdges);

  for (int i=1;i < 11;i++)
   {
    float em3=meanE_ch3_40[i];
    float em3b=meanE_ch3_40[8];
    mean40_3->SetBinContent(i+1,mean_ch3_40[i]);
    mean40_3->SetBinError(i+1,meanE_ch3_40[i]);
    TW3->SetBinContent(i+1,mean_ch3_40[8]-mean_ch3_40[i]);
    TW3->SetBinError(i+1,sqrt(em3*em3+em3b*em3b));
    rms40_3->SetBinContent(i+1,rms_ch3_40[i]);
    rms40_3->SetBinError(i+1,rmsE_ch3_40[i]);

    float em4=meanE_ch4_40[i];
    float em4b=meanE_ch4_40[8];
    mean40_4->SetBinContent(i+1,mean_ch4_40[i]);
    mean40_4->SetBinError(i+1,meanE_ch4_40[i]);
    TW4->SetBinContent(i+1,mean_ch4_40[8]-mean_ch4_40[i]);
    TW4->SetBinError(i+1,sqrt(em4*em4+em4b*em4b));
    rms40_4->SetBinContent(i+1,rms_ch4_40[i]);
    rms40_4->SetBinError(i+1,rmsE_ch4_40[i]);

    float em7=meanE_ch7_40[i];
    float em7b=meanE_ch7_40[8];
    mean40_7->SetBinContent(i+1,mean_ch7_40[i]);
    mean40_7->SetBinError(i+1,meanE_ch7_40[i]);
    TW7->SetBinContent(i+1,mean_ch7_40[8]-mean_ch7_40[i]);
    TW7->SetBinError(i+1,sqrt(em7*em7+em7b*em7b));
    rms40_7->SetBinContent(i+1,rms_ch7_40[i]);
    rms40_7->SetBinError(i+1,rmsE_ch7_40[i]);

    mean40_8->SetBinContent(i+1,mean_ch8_40[i]);
    mean40_8->SetBinError(i+1,meanE_ch8_40[i]);
    rms40_8->SetBinContent(i+1,rms_ch8_40[i]);
    rms40_8->SetBinError(i+1,rmsE_ch8_40[i]);
    float em8=meanE_ch8_40[i];
    float em8b=meanE_ch8_40[8];
    TW8->SetBinContent(i+1,mean_ch8_40[8]-mean_ch8_40[i]);
    TW8->SetBinError(i+1,sqrt(em8*em8+em8b*em8b));

    mean40_9->SetBinContent(i+1,mean_ch9_40[i]);
    mean40_9->SetBinError(i+1,meanE_ch9_40[i]);
    rms40_9->SetBinContent(i+1,rms_ch9_40[i]);
    rms40_9->SetBinError(i+1,rmsE_ch9_40[i]);
    float em9=meanE_ch9_40[i];
    float em9b=meanE_ch9_40[8];
    TW9->SetBinContent(i+1,mean_ch9_40[8]-mean_ch9_40[i]);
    TW9->SetBinError(i+1,sqrt(em9*em9+em9b*em9b));

   }
  mean40_9->SetLineColor(kYellow);
  rms40_9->SetLineColor(kYellow);
  TW9->SetLineColor(kYellow);

  mean40_8->SetLineColor(kRed);
  rms40_8->SetLineColor(kRed);
  TW8->SetLineColor(kRed);

  mean40_7->SetLineColor(kGreen);
  rms40_7->SetLineColor(kGreen);
  TW7->SetLineColor(kGreen);

  mean40_4->SetLineColor(kBlue);
  rms40_4->SetLineColor(kBlue);
  TW4->SetLineColor(kBlue);

  mean40_3->SetLineColor(kCyan);
  rms40_3->SetLineColor(kCyan);
  TW3->SetLineColor(kCyan);

  TCanvas c1;
  gStyle->SetOptStat(0);
  mean40_3->GetYaxis()->SetRangeUser(100.,700.);
  rms40_3->GetYaxis()->SetRangeUser(70.,230.);
  TW3->GetYaxis()->SetRangeUser(-40.,170.);
  mean40_3->Draw();
  mean40_4->Draw("same");
  mean40_7->Draw("same");
  mean40_8->Draw("same");
  mean40_9->Draw("same");
  auto leg=new TLegend(0.7,0.5,0.9,0.7);
  leg->SetHeader("Single pad ch.","C");
  leg->AddEntry(mean40_3,"channel 3","l");
  leg->AddEntry(mean40_4,"channel 4","l");
  leg->AddEntry(mean40_7,"channel 7","l");
  leg->AddEntry(mean40_8,"channel 8","l");
  leg->AddEntry(mean40_9,"channel 9","l");
  leg->Draw();
  c1.Print("plots/Run27-timeWalk-meanTime-vs-amplitude-singlePads-3-4-7-8-9.png");
  c1.Clear();
  rms40_3->Draw();
  rms40_4->Draw("same");
  rms40_7->Draw("same");
  rms40_8->Draw("same");
  rms40_9->Draw("same");
  auto leg2=new TLegend(0.7,0.7,0.9,0.9);
  leg2->SetHeader("Multipad ch.","C");
  leg2->AddEntry(rms40_3,"channel 3","l");
  leg2->AddEntry(rms40_4,"channel 4","l");
  leg2->AddEntry(rms40_7,"channel 7","l");
  leg2->AddEntry(rms40_8,"channel 8","l");
  leg2->AddEntry(rms40_9,"channel 9","l");
  leg2->Draw();
  c1.Print("plots/Run27-timeWalk-rmsTime-vs-amplitude-singlePads-3-4-7-8-9.png");
  c1.Clear();
  TW3->Draw();
  TW4->Draw("same");
  TW7->Draw("same");
  TW8->Draw("same");
  TW9->Draw("same");
  auto leg3=new TLegend(0.7,0.7,0.9,0.9);
  leg3->SetHeader("Multipad ch.","C");
  leg3->AddEntry(TW3,"channel 3","l");
  leg3->AddEntry(TW4,"channel 4","l");
  leg3->AddEntry(TW7,"channel 7","l");
  leg3->AddEntry(TW8,"channel 8","l");
  leg3->AddEntry(TW9,"channel 9","l");
  leg3->Draw();
  c1.Print("plots/Run27-timeWalk-correction-vs-amplitude-singlePads-3-4-7-8-9-45mV-ref.png");
  c1.Clear();
  gStyle->SetOptFit(1);
  auto walk=new TF1("walk","[0]/TMath::Sqrt(x+[1])+[2]",5.,150.);
  auto fex=new TF1("fex","expo(0)",5.,150.);
//  walk->SetParLimits(0,0.05,1000.);
//  walk->SetParLimits(1,-10.,10.);
//  walk->SetParLimits(2,-40.,40.);
  TW3->GetYaxis()->SetRangeUser(-20.,120.);
  TW3->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch3-45mV-ref.png");
  TW4->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch4-45mV-ref.png");
  TW7->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch7-45mV-ref.png");
  TW8->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch8-45mV-ref.png");
  TW9->Fit("fex");
  c1.Print("plots/Run27-timeWalk-correction-expfit5-ch9-45mV-ref.png");
  mean40_9->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-sing-ch9.png");
  mean40_8->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-sing-ch8.png");
  mean40_7->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-sing-ch7.png");
  mean40_4->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-sing-ch4.png");
  mean40_3->Fit("walk");
  c1.Print("plots/Run27-timeWalk-fit-sing-ch3.png");

 }
