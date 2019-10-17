#ifndef analysis_MyAnalysis_h
#define analysis_MyAnalysis_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "SampleAnalyzer/Interfaces/delphes/DelphesDataFormat.h"
#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"


namespace MA5
{
class MyAnalysis : public AnalyzerBase
{
  INIT_ANALYSIS(MyAnalysis,"MyAnalysis")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
   TH1F* myHisto1;
   TH1F* myHisto2;
};
}

#endif
