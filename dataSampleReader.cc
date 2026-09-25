#include <ROOT/RNTupleReader.hxx>
#include <TH1F.h>
#include <TFile.h>
#include <array>
#include <iostream>

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_171519_run44.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  auto hBase21=new TH1F("histBaseline21","Recalculated baseline, entries 1-20",200,-1.,1.);
  std::uint64_t numEntries = reader->GetNEntries();

  //event loop
  for (std::uint64_t i = 0; (i < numEntries)&& (i<30000); ++i) 
   {
    // viewArr(i) returns a const reference to the std::array<float, 64> at entry i
    const std::array<float, 64>& samples = viewArr(i);

    // Process elements inside the fixed-length vector
    float firstElement = samples[0];
    float lastElement = samples[63];
    float baseline=samples[1];
    for (int j=2;j<21;j++)
     baseline+=samples[j];
    baseline=baseline/20.;
    hBase21->Fill(baseline);
    if (!(i%1000)) 
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
  auto outFile=TFile::Open("baselines.root","RECREATE");
  hBase21->Write();
  outFile->Close();
 }
