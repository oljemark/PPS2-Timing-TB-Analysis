#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void percentileTuning()
 {
  auto h80=new TH1F("RMS_PMT_LGAD","Gaussian fit Sigma for interpolated 80% amplitude crossing time difference LGAD - PMT;SAMPIC channel;Sigma (ps)",10,-0.5,9.5);
  auto h70=new TH1F("RMS_PMT_LGAD7","Gaussian fit Sigma for interpolated 70% amplitude crossing time difference LGAD - PMT;SAMPIC channel;Sigma (ps)",10,-0.5,9.5);
  auto h50=new TH1F("RMS_PMT_LGAD8","Gaussian fit Sigma for interpolated 50% amplitude crossing time difference LGAD - PMT;SAMPIC channel;Sigma (ps)",10,-0.5,9.5);
  float sig50[10]={114.26, 113.792, 98.4176, 106.065, 105.76, 102.228, 103.388, 113.948, 116.729, 128.767};
  float sigE50[10]={1.63744, 0.855473, 0.639152, 1.11974, 1.11071, 0.458468, 0.383746, 1.47284, 1.41806, 1.90123};
  float sig70[10]={130.567, 119.127, 118.359, 110.276, 108.22, 112.793, 114.581, 104.375, 115.844, 120.785};
  float sigE70[10]={2.0752, 1.52378, 1.46971, 0.408467, 0.457176, 1.1648, 1.18417, 0.639238, 0.717245, 1.73962};
  float sig80[10]={118.184, 120.527, 109.881, 117.341 , 117.522, 113.071, 116.542, 119.772, 121.088, 129.404};
  float sigE80[10]={1.66237, 0.774273, 0.668781, 1.20062 , 1.21718, 0.471559, 0.437094, 1.49527, 1.58565, 2.00063};
  for (int i=0;i < 10;i++)
   {
    h50->SetBinContent(i+1,sig50[i]);
    h50->SetBinError(i+1,sigE50[i]);
    h70->SetBinContent(i+1,sig70[9-i]);
    h70->SetBinError(i+1,sigE70[9-i]);
    h80->SetBinContent(i+1,sig80[i]);
    h80->SetBinError(i+1,sigE80[i]);
   }
  h50->SetLineColor(kRed);
  h70->SetLineColor(kCyan);
  TCanvas c1;
  gStyle->SetOptStat(0);
  h50->Draw();
  h70->Draw("same");
  h80->Draw("same");
  auto leg=new TLegend(0.4,0.7,0.6,0.9);
  leg->AddEntry(h50,"half amplitude","l");
  leg->AddEntry(h70,"70% amplitude","l");
  leg->AddEntry(h80,"80% amplitude","l");
  leg->Draw();
  c1.Print("plots/Gauss-RMS-50-70-80thresh.png");
 }
