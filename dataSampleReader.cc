#include <ROOT/RNTupleReader.hxx>
#include <TH1F.h>
#include <TFile.h>
#include <array>
#include <iostream>

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_171519_run44.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  // may be int32 in 2025 data? No, in earlier version of SAMPIClyser before 0.2.0
  auto viewCh = reader->GetView<uint8_t>("Channel");
  auto hBase21=new TH1F("histBaseline21","Recalculated baseline, entries 1-20",2000,-1.,1.);
  auto hBase8=new TH1F("histBaseline21_ch8","Recalculated baseline, entries 1-20 (channel==8)",2000,-1.,1.);
  auto hAmpl8=new TH1F("histAmplitude21to53_ch8","Recalculated amplitude, entries 21-53 (channel==8)",2000,-1.,1.);
  auto hDC50ch8=new TH1F("histDC50sample_ch8","First sample reaching half amplitude, entries 21-53 (channel==8);Sample number (6400 MSamp/sec, 156 psec/sample)",64,-1.5,62.5);
  std::uint64_t numEntries = reader->GetNEntries();

  //event loop
  for (std::uint64_t i = 0; (i < numEntries)&& (i<130000); ++i)
   {
    // viewArr(i) returns a const reference to the std::array<float, 64> at entry i
    const std::array<float, 64>& samples = viewArr(i);
    const uint8_t ch=viewCh(i);

    // Process elements inside the fixed-length vector
    float firstElement = samples[0];
    float lastElement = samples[63];
    float baseline=samples[1];

    for (int j=2;j<21;j++)
     baseline+=samples[j];
    baseline=baseline/20.;
    hBase21->Fill(baseline);
    if (ch==8)
     {
      hBase8->Fill(baseline);
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
      hDC50ch8->Fill(j50);
      hAmpl8->Fill(amplitude);
     }
    if (!(i%1000)) 
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
  auto outFile=TFile::Open("amplitudes.root","RECREATE");
  hBase21->Write();
  hBase8->Write();
  hAmpl8->Write();
  hDC50ch8->Write();
  outFile->Close();
 }
