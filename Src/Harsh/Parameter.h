#pragma once

#include <vector>
#include <string>

using std::vector;
using std::string;

class Parameter {

public:
  Parameter();

  void setVerbose(uint8_t verbose);
  void setMaxIter(int32_t maxIter);
  void setSubIter(int32_t subIter);
  void setChunkSize(int32_t chunkSize);
  void setOverlapSize(int32_t overlapSize);
  void setSampleSize(int32_t sampleSize);

  void setSeqFile(const char* sFile);
  void setRefFile(const char* rFile);
  void setMapFile(const char* mFile);
  void setOutFile(const char* oFile);
 
  void setSeqError(double seqError);
  void setEmitError(double emitError);
  void setTranNe(double transNe);
  void setSmoothMu(double smoothMu);

  uint8_t getVerbose() const;
  int32_t getMaxIter() const;
  uint32_t getSubIter() const;
  bool getHasRef() const;
  int32_t getChunkSize() const;
  int32_t getOverlapSize() const;  
  int32_t getSampleSize() const;  

  const char* getSeqFile() const;
  const char* getRefFile() const;
  const char* getMapFile() const;
  const char* getOutFile() const;
  
  double getSeqError() const;
  double getSeqRight() const;
  double getSeqLogError() const;
  double getSeqLogRight() const;
  double getEmitError() const;
  double getEmitRight() const;
  double getEmitLogError() const;
  double getEmitLogRight() const;

  double getTranNe() const;
  double getTranLogDiff(uint32_t index) const;
  double getTranLogSame(uint32_t index) const;
  double getTranLogPlus(uint32_t index) const;
  double getSmoothMu() const;

  double getForwardRecombProb(uint32_t snpIndex) const;
  double getBackwardRecombProb(uint32_t snpIndex) const;
  char getMajorAllele(uint32_t snpIndex) const;
  char getMinorAllele(uint32_t snpIndex) const;
  
  uint32_t getNSnp() const;
  void initFromMapFile();
  int validateLine(const string& line) const;
  void renewTran(uint32_t nHap);

private:
  
  /* Program Parameters */
  uint8_t verbose_;
  int32_t maxIter_;
  uint32_t subIter_;
  bool hasRef_;
  int32_t chunkSize_;
  int32_t overlapSize_;
  int32_t sampleSize_;

  /* Files */
  const char* sFile_;
  const char* rFile_;
  const char* mFile_;
  const char* oFile_;

  /* Model Parameters */
  double seqError_;
  double emitError_;
  double tranNe_;
  double smoothMu_;
  uint32_t nHap_;

  /* From SNP Map File*/
  uint32_t nSnp_;
  vector<string> chrom_;
  vector<string> rs_;
  vector<double> gmap_;
  vector<uint32_t> map_;
  vector<char> major_;
  vector<char> minor_;

  /* precomputed parameters */
  double logSeqError_;
  double logSeqRight_;
  double logEmitError_;
  double logEmitRight_;
  vector<double> tranLogDiff_;
  vector<double> tranLogSame_;
  vector<double> tranLogPlus_;
  
};
