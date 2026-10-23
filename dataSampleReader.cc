#include <ROOT/RNTupleReader.hxx>
#include <TH1F.h>
#include <TH1D.h>
#include <TFile.h>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

struct waves
 {
  int64_t index{0}; // event number
  int64_t coarseTime{0}; // time to first sample, ps
  float amplitude{0.};
  int64_t timeAtDC50{0}; // time at midpoint of rising edge, ps (tune to 40%)
  int maxbin{0}; // max amplitude bin, out of 64
  int binDC50{0}; // 1st half amplitude bin, out of 64
  float finebinDC50{0.}; // interpolate half amplitude bin position, between binDC50 and binDC50-1
 };

void dataSampleReader()
 {
  auto reader = ROOT::RNTupleReader::Open("sampic_hits", "root_20260810_run44_full.root");
  auto viewArr = reader->GetView<std::array<float, 64>>("DataSample");
  // may be int32 in 2025 data? No, in earlier version of SAMPIClyser before 0.2.0
  auto viewCh = reader->GetView<uint8_t>("Channel");
  auto viewTSPSint = reader->GetView<int64_t>("FirstSampleTime_in_ps");
  auto hBase21=new TH1F("histBaseline21","Recalculated baseline, entries 1-20",2000,-1.,1.);
  auto hPolar=new TH1F("histPolarity","Recalculated amplitude per SAMPIC channel (negative numbers are channels with amplitude < 0);SAMPIC channel",65,-32.5,32.5);

  auto hBase32=new TH1F("histBaselinePMT","Recalculated baseline, entries 1-20 (channel 32=MCP-PMT)",2000,-1.,1.);
  auto hBaseRms32=new TH1F("histBaselineRmsPMTTail6sig","Recalculated baseline rms, entries 1-20 (channel 32=MCP-PMT);Baseline rms above 6-sigma cut (mV)",2000,1.,-1.);
  auto hAmpl32=new TH1F("histAmplitudePMT","Recalculated amplitude (negative), entries 21-53 (channel 32=MCP-PMT)",2000,-1.,1.);
  auto hDC50_32=new TH1F("histDC40SamplePMT","First sample to reach 40 percent of amplitude, entries 21-53 (channel 32=MCP-PMT);Sample number",64,-1.5,62.5);

  auto hNearC=new TH1F("histNearCh_coarse","Hits in channels 0-9, with firstSampleTime within 10 ns",1,0.,1.);
  auto hNearFine=new TH1F("histNearCh_fine","Hits in channels 0-9, with rising edge DC50Time within 5 ns",1,0.,1.);
  auto hNearPMT=new TH1F("histNearPMT_fine","Hits in channels 0-9, with rising edge DC50Time within 1400 ps of PMT-MCP channel 32",1,0.,1.);
  int64_t lastDC50[10]={0}; // latest picosec time (DC50 corrected) for channel 0-9
  std::vector<TH1F*> hBases,hBaseRms,hAmpls,hAmp2,hDC50s,hCombos,hComboPMT,hWalkComboPMT,hBipole,hPosPeaks;
  std::vector<TH1F*> hWavef8max,hWavef8min;
  std::vector<TH1F*> hBigAmp,h2ndAmp;
  std::vector<TH1D*> hRealTrkIndexDiffPMT;
  hRealTrkIndexDiffPMT.reserve(10);
  hBigAmp.reserve(10);
  h2ndAmp.reserve(10);
  hBases.reserve(10);
  hBaseRms.reserve(10);
  hAmpls.reserve(10);
  hPosPeaks.reserve(10);
  hBipole.reserve(11);
  hAmp2.reserve(10);
  hDC50s.reserve(10);
  hComboPMT.reserve(10);
  hWalkComboPMT.reserve(100);
  hCombos.reserve(100);
  hWavef8max.reserve(200);
  hWavef8min.reserve(200);
  std::vector<waves> chWaves[10];
  std::vector<waves> pmtWaves;
  std::map<uint8_t,int> noisyBase;
  std::map<uint8_t,float> noiseThr;
  noiseThr[0]=1.97;
  noiseThr[1]=2.01;
  noiseThr[2]=1.99;
  noiseThr[3]=1.86;
  noiseThr[4]=2.06;
  noiseThr[5]=1.97;
  noiseThr[6]=1.95;
  noiseThr[7]=1.81;
  noiseThr[8]=1.73;
  noiseThr[9]=1.83;
  noiseThr[32]=3.43; // fitted gaussian in the tail above 1.65 (6 sigma from central peak)

  noisyBase[32]=0;
  for (uint8_t ch=0;ch<10;ch++)
   noisyBase[ch]=0;
  for (uint16_t ch=0; ch<200;ch++)
  {
     std::string name="waveform_ch_neg_"+std::to_string(ch);
     std::string title="Baseline high rms (>2 mV) example waveform in channels 0-9;Sample number;Amplitude above baseline (mV)";
     hWavef8min.push_back(new TH1F(name.c_str(),title.c_str(),65,-0.5,64.5));
     name="waveform_ch_pos_"+std::to_string(ch);
     title="Positive amplitude above baseline RMS noise threshold (6 sigma), example waveform in channels 0-9;Sample number;Amplitude above baseline (mV)";
     hWavef8max.push_back(new TH1F(name.c_str(),title.c_str(),65,-0.5,64.5));
  }
   for (uint8_t ch=0; ch<11;ch++)
    {
     std::string name="smallerPolarityAmplitude_"+std::to_string(ch);
     if (ch==10)
      name="smallerPolarityAmplitude_32";
     std::string title="Smaller of amplitude above or below baseline for LGAD ch"+ std::to_string(ch)+";Smaller amplitude (mV)";
     if (ch==10)
      title="Smaller of amplitude above or below baseline for PMT ch 32;Smaller amplitude (mV)";
     hBipole.push_back(new TH1F(name.c_str(),title.c_str(),500,1.,-1.));
    }

   for (uint8_t ch2=0; ch2<10;ch2++)
    {
     const std::string amplis[11]={"2-5mV","5-10mV","10-15mV","15-20mV","20-25mV","25-30mV","30-35mV","35-45mV","45-60mV","60-80mV","80-150mV"};
     std::string name="timeDiffPMT_ch32_"+std::to_string(ch2);
     std::string title="Time difference of first sample reaching 40 pct of (negative) amplitude (max 1400ps),  PMT channel vs LGAD ch"+ std::to_string(ch2)+";Time difference (ps)";
     if (ch2==8)
      {
       hComboPMT.push_back(new TH1F(name.c_str(),title.c_str(),320,-200.,1400.));
      }
     else
     {
      hComboPMT.push_back(new TH1F(name.c_str(),title.c_str(),320,-200.,1400.));
     }
    for (uint8_t amp=0; amp<11;amp++)
     {
      name="timeDiffPMT_ch32_"+std::to_string(ch2)+"_"+amplis[amp];
      title="Time difference of first sample reaching 40 pct of amplitude (max 1400ps),  PMT channel vs LGAD ch"+ std::to_string(ch2)+", max amplitude in range "+amplis[amp]+";Time difference (ps)";
      hWalkComboPMT.push_back(new TH1F(name.c_str(),title.c_str(),80,-200.,1400.));
     }

     }
  for (uint8_t ch=0; ch<10;ch++)
   {
    const float binEdges[12]={2,5,10,15,20,25,30,35,45,60,80,150};

    std::string name="histBaseline_ch" + std::to_string(ch);
    std::string title="Recalculated baseline, entries 1-20, channel=" + std::to_string(ch);
    hBases.push_back(new TH1F(name.c_str(),title.c_str(),2000,-1.,1.));

    name="histPMTIndexOverlap_ch" + std::to_string(ch);
    title="Event number difference between LGAD and coincidence in PMT (ev_LGAD - ev_PMT), channel=" + std::to_string(ch);
    hRealTrkIndexDiffPMT.push_back(new TH1D(name.c_str(),title.c_str(),2000,1.,0.));

    name="histPosPeaks_ch" + std::to_string(ch);
    title="Count of full positive peaks (rising+max+falling edges), channel=" + std::to_string(ch);
    hPosPeaks.push_back(new TH1F(name.c_str(),title.c_str(),40,-0.5,39.5));

    name="histBaselineVar_ch" + std::to_string(ch);
    title="RMS of recalculated baseline, entries 1-20, channel=" + std::to_string(ch)+";Baseline standard deviation (mV)";
    hBaseRms.push_back(new TH1F(name.c_str(),title.c_str(),2000,1.,-1.));

    name="histBiggestPeakAmp_ch" + std::to_string(ch);
    title="Biggest peak amplitude above recalculated baseline, entries 1-20, channel=" + std::to_string(ch)+";Amplitude over waveform baseline (mV)";
    hBigAmp.push_back(new TH1F(name.c_str(),title.c_str(),2000,1.,-1.));

    name="hist2ndBiggestPeakAmp_ch" + std::to_string(ch);
    title="Second biggest peak amplitude above recalculated baseline (may be negative), entries 1-20, channel=" + std::to_string(ch)+";Amplitude over waveform baseline (mV)";
    h2ndAmp.push_back(new TH1F(name.c_str(),title.c_str(),2000,1.,-1.));

    name="histAmplitude21to53_ch"+ std::to_string(ch);
    title="Recalculated amplitude, entries 21-53, channel="+ std::to_string(ch);
    hAmpls.push_back(new TH1F(name.c_str(),title.c_str(),2000,-1.,1.));
    name=name+"varBinning";
    title="Recalculated amplitude (mV), entries 21-53, channel="+ std::to_string(ch)+";Amplitude (mV)";
    hAmp2.push_back(new TH1F(name.c_str(),title.c_str(),11,binEdges));

    name="histDC40sample_ch"+ std::to_string(ch);
    title="First sample reaching 40 percent of amplitude, entries 21-53, channel="+ std::to_string(ch)+";Sample number";
    hDC50s.push_back(new TH1F(name.c_str(),title.c_str(),64,-1.5,62.5));

    for (uint8_t ch2=0; ch2<10;ch2++)
     {
      name="timeDiff2_ch"+ std::to_string(ch)+"_"+std::to_string(ch2);
      title="Time difference of first sample reaching 40 pct of amplitude (max 5ns),  channel="+ std::to_string(ch)
           +" vs "+ std::to_string(ch2)+";Time difference (ps)";
      hCombos.push_back(new TH1F(name.c_str(),title.c_str(),400,-5000.,5000.));

     }
   }
  auto hPs50ch8=new TH1D("histPicosec40sample_ch8","First sample time in integer picoseconds reaching 40 pct of amplitude, entries 21-53 (channel==8);Time (ps)",200,0.,8e15);
  std::uint64_t numEntries = reader->GetNEntries();
  int neg8=0;
  int pos8=0;

  //event loop
  for (std::uint64_t i = 0; (i < numEntries)&& (i<13000000); ++i)
   {
    // viewArr(i) returns a const reference to the std::array<float, 64> at entry i
    const std::array<float, 64>& samples = viewArr(i);
    const uint8_t ch=viewCh(i);
    const float thr=0.4;
    const int64_t picosec=viewTSPSint(i); // first sample time, rounded to nearest ps

    // Process elements inside the fixed-length vector
    float firstElement = samples[0];
    float lastElement = samples[63];
    float baseline=samples[1];

    for (int j=2;j<21;j++)
     baseline+=samples[j];
    baseline=baseline/20.;
    hBase21->Fill(baseline);
    float var=0.;
    for (int j=1;j<21;j++)
     var+=(samples[j]-baseline)*(samples[j]-baseline);
    var=var/20.;
    const float rms=std::sqrt(var);
    if (1000.*rms>noiseThr[ch])
     noisyBase[ch]++;

    if (ch==32)
     {
      waves pmt;
      pmt.index=i;
      pmt.coarseTime = picosec;
      hBase32->Fill(baseline);
      if (1000.*rms > noiseThr[32])
       hBaseRms32->Fill(1000.*rms);
      float min=baseline;
      float max=baseline;
      int jmin=20;
      for (int j=21;j<54;j++)
       {
        if (max<samples[j])
         max=samples[j];
        if (min>samples[j])
         {
          min=samples[j];
          jmin=j;
         }
       }
      float amplitude=min-baseline;
      float PAmplitude=max-baseline;
      const int pole=((PAmplitude > -amplitude) ? 32 : -32);
      hPolar->Fill(pole);
      const float lesser = (pole<0 ? PAmplitude : -amplitude);
      hBipole[10]->Fill(1000.*lesser);
      pmt.amplitude=amplitude;
      hAmpl32->Fill(amplitude);
      pmt.maxbin=jmin;
      float dc50=amplitude*thr;
      int j50=20;
      while (j50<jmin)
       {
        if (samples[j50]<baseline+dc50)
         break;
        j50++;
       }
      float upper50=(samples[j50]-baseline)/amplitude;
      float lower50=(samples[j50-1]-baseline)/amplitude;
      if (amplitude < 0)
       {
        float subrange=upper50-lower50;
        float downer=upper50 - thr;
        float delta=downer/subrange;
        pmt.finebinDC50=j50 - delta;
       }
      else
      {
       pmt.finebinDC50=j50+0.;
      }
      pmt.binDC50=j50;
      hDC50_32->Fill(j50);
      const int ps50=(pmt.finebinDC50 * 1000) / 6.4;
      const int64_t time50ps = picosec + ps50;
      pmt.timeAtDC50=time50ps;
      pmtWaves.push_back(pmt);
     }
    if (ch<10)
     {
      waves wave;
      wave.index=i;
      wave.coarseTime = picosec;
      hBases[ch]->Fill(baseline);
      hBaseRms[ch]->Fill(1000.*rms);
      float max=baseline;
      float min=baseline;
      int rising[64]={0};
      int falling[64]={0};
      int locmin[64]={0};
      int locmax[64]={0};
      int posPeaks=0;
      std::vector<float> amps;
      enum states {kNA,kRise,kFall,kMin,kMax};
      int jmax=20;
      for (int j=21;j<54;j++)
       {
        if (min>samples[j])
         min=samples[j];
        if (max<samples[j])
         {
          max=samples[j];
          jmax=j;
         }
       }
      float amplitude=max-baseline;
      float MAmplitude=min-baseline;
      const int pole=((-MAmplitude > amplitude) ? -ch : ch);
      hPolar->Fill(pole);
      const float lesser = (pole>0 ? -MAmplitude : amplitude);
      hBipole[ch]->Fill(1000.*lesser);

      for (int j=21;j<54;j++)
       {
        float ldiff1=samples[j-2]-samples[j-1];
        float ldiff2=samples[j-1]-samples[j];
        float rdiff1=samples[j+2]-samples[j+1];
        float rdiff2=samples[j+1]-samples[j];
        float meanleft=(ldiff1+ldiff2)/2.;
        float meanright=(rdiff1+rdiff2)/2.;
        if ((meanleft<0)&&(meanright>0))
         rising[j]=1;
        else
         rising[j]=2;
        if ((meanleft>0)&&(meanright<0))
         falling[j]=1;
        else
         falling[j]=2;
        if ((meanleft>0)&&(meanright>0))
         locmin[j]=1;
        else
         locmin[j]=2;
        if ((meanleft<0)&&(meanright<0))
         locmax[j]=1;
        else
         locmax[j]=2;
       }
      states st=kNA;
      float amp=-99.;
      for (int j=21;j<54;j++)
       {
        if (rising[j]==1)
         st=kRise;
        if ((locmax[j]==1)&&(st==kRise))
         {
          st=kMax;
          amp=1000.*(samples[j]-baseline);
         }
        if ((falling[j]==1)&&(st==kMax))
         {
          st=kFall;
          posPeaks++;
          amps.push_back(amp);
         }
       }
      hPosPeaks[ch]->Fill(posPeaks+0.);
      if (!amps.empty())
       {
        std::sort(amps.begin(),amps.end(),std::greater<float>());
        if (posPeaks>0)
         {
          hBigAmp[ch]->Fill(amps.at(0));
         }
        if ((posPeaks>1) && (amps.size()>1))
         h2ndAmp[ch]->Fill(amps.at(1));
       }
      if ((neg8<200)&&(rms*1000.>noiseThr[ch]))
      {
       for (int j=1;j<54;j++)
         hWavef8min[neg8]->SetBinContent(j+1,1000.*(samples[j]-baseline));
       neg8++;
      }
      if ((pos8<200)&&(amplitude*1000.>noiseThr[ch]))
      {
       for (int j=1;j<54;j++)
         hWavef8max[pos8]->SetBinContent(j+1,1000.*(samples[j]-baseline));
       pos8++;
      }

      wave.amplitude=amplitude;
      wave.maxbin=jmax;
      float dc50=amplitude*thr;
      int j50=20;
      while (j50<jmax)
       {
        if (samples[j50]>baseline+dc50)
         break;
        j50++;
       }
      float upper50=(samples[j50]-baseline)/amplitude;
      float lower50=(samples[j50-1]-baseline)/amplitude;
      if (amplitude > 0)
       {
        float subrange=upper50-lower50;
        float downer=upper50 - thr;
        float delta=downer/subrange;
        wave.finebinDC50=j50 - delta;
       }
      else
      {
       wave.finebinDC50=j50+0.;
      }
      wave.binDC50=j50;
      const int ps50=(wave.finebinDC50 * 1000) / 6.4;
      const int64_t time50ps = picosec + ps50;
      wave.timeAtDC50=time50ps;
      chWaves[ch].push_back(wave);
      if ((amplitude > 0) && (jmax > 20))
       {
        // for (uint8_t ch2=0; ch2<10; ch2++)
        //  hCombos[10*ch+ch2]->Fill(time50ps - lastDC50[ch2]);
        lastDC50[ch]=time50ps;
       }
      hDC50s[ch]->Fill(j50);
      if (ch == 8)
       hPs50ch8->Fill(time50ps);
      if (!(i%100))
       std::cout<<"time at 40 pct amplitude(ps):"<<time50ps<<"ps j50:"<<j50<<" ch:"<<std::to_string(ch)<<std::endl;
      hAmpls[ch]->Fill(amplitude);
      hAmp2[ch]->Fill(1000.*amplitude);
     }
    if (!(i%10000))
     {
      std::cout<<"i"<<i<<"First entry: "<<firstElement<<" last:"
      <<lastElement<<std::endl;
     }
   }
  for (const auto& [cha, entries] : noisyBase)
   std::cout<<"Channel "+std::to_string(cha)+", entries:"<<entries<<std::endl;
  auto outFile=TFile::Open("maps-run044-40pct-v1.root","RECREATE");
  hBase21->Write();
  hBase32->Write();
  hBaseRms32->Write();
  hAmpl32->Write();
  hDC50_32->Write();
  std::cout<<"MCP-PMT triggers:"<<pmtWaves.size()<<std::endl;
  for (uint16_t ch=0; ch<11;ch++)
   hBipole[ch]->Write();
  for (uint16_t ch8=0; ch8<200;ch8++)
  {
   hWavef8max[ch8]->Write();
   hWavef8min[ch8]->Write();
  }
  for (uint8_t ch=0; ch<10;ch++)
   {
    std::cout<<chWaves[ch].size()<<std::endl;
    hBases[ch]->Write();
    hBigAmp[ch]->Write();
    h2ndAmp[ch]->Write();
    hPosPeaks[ch]->Write();
    hBaseRms[ch]->Write();
    hAmpls[ch]->Write();
    hAmp2[ch]->Write();
    hDC50s[ch]->Write();
    std::string name="data/coincidencePerChannel"+std::to_string(ch)+".txt";
    std::ofstream chou(name.c_str(),std::ofstream::trunc);
    bool opn=chou.is_open();
    for (const auto& pmt: pmtWaves)
     {
      for (const auto& wave1: chWaves[ch])
       {
        if (std::abs(wave1.timeAtDC50 - pmt.timeAtDC50) < (ch==8 ? 1400 : 1400))
         {
          std::string labl=std::to_string(ch)+",32";
          hNearPMT->Fill(labl.c_str(),1.);
          hRealTrkIndexDiffPMT[ch]->Fill(wave1.index-pmt.index);
          if (opn)
           chou<<(wave1.index-pmt.index)<<" ("<<wave1.index<<" "<<pmt.index<<")"<<std::endl;
          hComboPMT[ch]->Fill(wave1.timeAtDC50 - pmt.timeAtDC50);
          const auto mV=1000.*wave1.amplitude;
          if ((mV>2.)&&(mV<150.))
           {
            uint8_t amp=0;
            if (mV<5.)
             amp=0;
            else
             {
              if (mV<35.)
               amp=1+(mV - 5.)/5;
              else
               {
                if (mV<45.)
                 amp=7;
                else
                 {
                  if (mV<60.)
                   amp=8;
                  else
                   {
                    if (mV<80.)
                     amp=9;
                    else
                     {
                      if (mV<150.)
                       amp=10;
                     }
                  }
               }
             }
           }
            hWalkComboPMT[11*ch+amp]->Fill(wave1.timeAtDC50 - pmt.timeAtDC50);
           }
         }
       }
     }
    if (opn)
     chou.close();
    hRealTrkIndexDiffPMT[ch]->Write();
    hComboPMT[ch]->Write();
    for (uint8_t amp=0; amp<11;amp++)
      hWalkComboPMT[ch*11+amp]->Write();

    for (uint8_t ch2=0; ch2<10;ch2++)
     {
      if (ch2==ch) continue;
      for (const auto& wave1: chWaves[ch])
       {
        for (const auto& wave2: chWaves[ch2])
         {
          if (std::abs(wave1.timeAtDC50 - wave2.timeAtDC50) < 5000)
           {
            std::string labl=std::to_string(ch)+","+std::to_string(ch2);
            hNearFine->Fill(labl.c_str(),1.);
            hCombos[10*ch+ch2]->Fill(wave1.timeAtDC50 - wave2.timeAtDC50);
           }
         }
       }
      hCombos[10*ch+ch2]->Write();
     }
   }
  hPolar->Write();
  hNearC->Write();
  hNearFine->Write();
  hNearPMT->Write();
  hPs50ch8->Write();
  outFile->Close();
 }
