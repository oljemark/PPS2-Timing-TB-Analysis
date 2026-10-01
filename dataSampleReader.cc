#include <ROOT/RNTupleReader.hxx>
#include <TH1F.h>
#include <TH1D.h>
#include <TFile.h>
#include <array>
#include <vector>
#include <string>
#include <iostream>

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_run44_full.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  // may be int32 in 2025 data? No, in earlier version of SAMPIClyser before 0.2.0
  auto viewCh = reader->GetView<uint8_t>("Channel");
  auto viewTSPSint = reader->GetView<int64_t>("FirstSampleTime_in_ps");
  auto hBase21=new TH1F("histBaseline21","Recalculated baseline, entries 1-20",2000,-1.,1.);
  std::vector<TH1F*> hBases,hAmpls,hDC50s;
  hBases.reserve(10);
  hAmpls.reserve(10);
  hDC50s.reserve(10);
  for (uint8_t ch=0; ch<10;ch++)
   {
    std::string name="histBaseline_ch" + std::to_string(ch);
    std::string title="Recalculated baseline, entries 1-20, channel=" + std::to_string(ch);
    hBases.push_back(new TH1F(name.c_str(),title.c_str(),2000,-1.,1.));

    name="histAmplitude21to53_ch"+ std::to_string(ch);
    title="Recalculated amplitude, entries 21-53, channel="+ std::to_string(ch);
    hAmpls.push_back(new TH1F(name.c_str(),title.c_str(),2000,-1.,1.));

    name="histDC50sample_ch"+ std::to_string(ch);
    title="First sample reaching half amplitude, entries 21-53, channel="+ std::to_string(ch)+";Sample number";
    hDC50s.push_back(new TH1F(name.c_str(),title.c_str(),64,-1.5,62.5));
   }
  auto hPs50ch8=new TH1D("histPicosec50sample_ch8","First sample time in integer picoseconds reaching half amplitude, entries 21-53 (channel==8);Time (ps)",200,0.,8e15);
  std::uint64_t numEntries = reader->GetNEntries();

  //event loop
  for (std::uint64_t i = 0; (i < numEntries)&& (i<13000000); ++i)
   {
    // viewArr(i) returns a const reference to the std::array<float, 64> at entry i
    const std::array<float, 64>& samples = viewArr(i);
    const uint8_t ch=viewCh(i);
    const int64_t picosec=viewTSPSint(i); // first sample time, rounded to nearest ps

    // Process elements inside the fixed-length vector
    float firstElement = samples[0];
    float lastElement = samples[63];
    float baseline=samples[1];

    for (int j=2;j<21;j++)
     baseline+=samples[j];
    baseline=baseline/20.;
    hBase21->Fill(baseline);
    if (ch<10)
     {
      hBases[ch]->Fill(baseline);
      float max=baseline;
      int jmax=20;
      for (int j=21;j<54;j++)
       if (max<samples[j])
        {
         max=samples[j];
         jmax=j;
        }
      float amplitude=max-baseline;
      float dc50=amplitude*0.5;
      int j50=20;
      while (j50<jmax)
       {
        if (samples[j50]>baseline+dc50)
         break;
        j50++;
       }
      const int ps50=(j50 * 1000) / 6.4;
      const int64_t time50ps = picosec + ps50;
      hDC50s[ch]->Fill(j50);
      if (ch == 8)
       hPs50ch8->Fill(time50ps);
      if (!(i%100))
       std::cout<<"time at half amplitude(ps):"<<time50ps<<"ps j50:"<<j50<<" ch:"<<std::to_string(ch)<<std::endl;
      hAmpls[ch]->Fill(amplitude);
     }
    if (!(i%10000))
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
  auto outFile=TFile::Open("amplitudes-run044-tenCh.root","RECREATE");
  hBase21->Write();
  for (uint8_t ch=0; ch<10;ch++)
   {
    hBases[ch]->Write();
    hAmpls[ch]->Write();
    hDC50s[ch]->Write();
   }
  hPs50ch8->Write();
  outFile->Close();
 }
