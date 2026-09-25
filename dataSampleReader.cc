#include <ROOT/RNTupleReader.hxx>
#include <TH1F.h>
#include <TFile.h>
#include <array>
#include <iostream>

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_171519_run44.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  // may be int32 in 2025 data?
  auto viewCh = reader->GetView<uint8_t>("Channel");
  auto hBase21=new TH1F("histBaseline21","Recalculated baseline, entries 1-20",2000,-1.,1.);
  auto hBase8=new TH1F("histBaseline21_ch8","Recalculated baseline, entries 1-20 (channel==8)",2000,-1.,1.);
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
     hBase8->Fill(baseline);
    if (!(i%1000)) 
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
  auto outFile=TFile::Open("baselines.root","RECREATE");
  hBase21->Write();
  hBase8->Write();
  outFile->Close();
 }
