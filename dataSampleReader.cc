#include <ROOT/RNTupleReader.hxx>
#include <array>
#include <iostream>

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_171519_run44.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  std::uint64_t numEntries = reader->GetNEntries();

  //event loop
  for (std::uint64_t i = 0; (i < numEntries)&& (i<30000); ++i) 
   {
    // viewArr(i) returns a const reference to the std::array<float, 64> at entry i
    const std::array<float, 64>& currentArray = viewArr(i);

    // Process elements inside the fixed-length vector
    float firstElement = currentArray[0];
    float lastElement = currentArray[63];
    if (!(i%1000)) 
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
 }
