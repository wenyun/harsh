#include "Parameter.h"

#include "HarshIO.h"

#include <fstream>
#include <sstream>
#include <cmath>

using std::ifstream;
using std::istringstream;

extern char line[1000];

Parameter::Parameter() : 
  verbose_(1),
  maxIter_(1000),
  subIter_(100),
  hasRef_(false),
  chunkSize_(5000),
  overlapSize_(500),
  sampleSize_(30),
  sFile_(""),
  rFile_(""),
  mFile_(""),
  oFile_(""),
  seqError_(0.01),
  emitError_(2e-3),
  tranNe_(0.001),
  smoothMu_(1.0),
  nHap_(0),
  nSnp_(0),
  chrom_(0),
  rs_(0),
  gmap_(0),
  map_(0),
  major_(0),
  minor_(0) {

  logSeqError_ = log(seqError_);
  logSeqRight_ = log(1 - seqError_);
  logEmitError_ = log(emitError_);
  logEmitRight_ = log(1 - emitError_);
}



void Parameter::setVerbose(uint8_t verbose) {
  verbose_ = verbose;
}

void Parameter::setMaxIter(int32_t maxIter) {
  maxIter_ = maxIter;
}

void Parameter::setSubIter(int32_t subIter) {
  subIter_ = subIter;
}

void Parameter::setChunkSize(int32_t chunkSize) {
  chunkSize_ = chunkSize;
}

void Parameter::setOverlapSize(int32_t overlapSize) {
  overlapSize_ = overlapSize;
}

void Parameter::setSampleSize(int32_t sampleSize) {
  sampleSize_ = sampleSize;
}

void Parameter::setSeqError(double seqError) {
  seqError_ = seqError;
  logSeqError_ = log(seqError_);
  logSeqRight_ = log(1 - seqError_);
}

void Parameter::setEmitError(double emitError) {
  emitError_ = emitError;
  logEmitError_ = log(emitError_);
  logEmitRight_ = log(1 - emitError_);
}

void Parameter::setTranNe(double tranNe) {
  tranNe_ = tranNe;
}

void Parameter::setSmoothMu(double smoothMu) {
  smoothMu_ = smoothMu;
}

void Parameter::setSeqFile(const char* sFile) {
  ifstream fin(sFile);
  if (fin.good()) {
    sFile_ = sFile;
  } else {
    sprintf(line,
            "Sequencing file %s does not exist\n",
            sFile);
    harshErrorExit(line);
  }
}

void Parameter::setRefFile(const char* rFile) {
  ifstream fin(rFile);
  if (fin.good()) {
    rFile_ = rFile;
    hasRef_ = true;
  } else {
    sprintf(line,
            "Reference file %s does not exist, swithcing to haplotype assembly\n",
            rFile);
    harshWarning(line);
    hasRef_ = false;
  }
}

void Parameter::setMapFile(const char* mFile) {
  ifstream fin(mFile);
  if (fin.good()) {
    mFile_ = mFile;
  } else {
    sprintf(line,
            "SNP map file %s does not exist\n",
            mFile);
    harshErrorExit(line);
  }
}

void Parameter::setOutFile(const char* oFile) {
    oFile_ = oFile;
}

uint8_t Parameter::getVerbose() const {
  return verbose_;
}

int32_t Parameter::getMaxIter() const {
  return maxIter_;
}

uint32_t Parameter::getSubIter() const {
  return subIter_;
}

bool Parameter::getHasRef() const {
  return hasRef_;
}

int32_t Parameter::getChunkSize() const {
  return chunkSize_;
}

int32_t Parameter::getOverlapSize() const {
  return overlapSize_;
}

int32_t Parameter::getSampleSize() const {
  return sampleSize_;
}

double Parameter::getSeqError() const {
  return seqError_;
}

double Parameter::getSeqRight() const {
  return 1.0 - seqError_;
}

double Parameter::getSeqLogError() const {
  return logSeqError_;
}

double Parameter::getSeqLogRight() const {
  return logSeqRight_;
}

double Parameter::getEmitError() const {
  return emitError_;
}

double Parameter::getEmitRight() const {
  return 1.0 - emitError_;
}

double Parameter::getEmitLogError() const {
  return logEmitError_;
}

double Parameter::getEmitLogRight() const {
  return logEmitRight_;
}

double Parameter::getTranNe() const {
  return tranNe_;
}

double Parameter::getTranLogDiff(uint32_t index) const {
  return tranLogDiff_[index];
}

double Parameter::getTranLogSame(uint32_t index) const {
  return tranLogSame_[index];
}

double Parameter::getTranLogPlus(uint32_t index) const {
  return tranLogPlus_[index];
}

double Parameter::getSmoothMu() const {
  return smoothMu_;
}

const char* Parameter::getSeqFile() const {
  return sFile_;
}

const char* Parameter::getRefFile() const {
  return rFile_;
}

const char* Parameter::getMapFile() const {
  return mFile_;
}

const char* Parameter::getOutFile() const {
  return oFile_;
}

char Parameter::getMajorAllele(uint32_t snpIndex) const {
  if (snpIndex < nSnp_) {
    return major_[snpIndex];
  } else {
    return '~';
  }
}

char Parameter::getMinorAllele(uint32_t snpIndex) const {
  if (snpIndex < nSnp_) {
    return minor_[snpIndex];
  } else {
    return '~';
  }
}

uint32_t Parameter::getNSnp() const {
  return nSnp_;
}

void Parameter::initFromMapFile() {
  
  ifstream fin(mFile_);
  string line;
  string chrom;
  string rs;
  double gmap;
  uint32_t map;
  char minor;
  char major;

  while(getline(fin, line)) {
    
    int nToken = validateLine(line);
    if (nToken == 6) {
      istringstream sin(line);
      sin >> chrom >> rs >> gmap >> map >> major >> minor;
      
      chrom_.push_back(chrom);
      rs_.push_back(rs);
      gmap_.push_back(gmap);
      map_.push_back(map);
      major_.push_back(major);
      minor_.push_back(minor);
    } else {
      harshErrorExit("Map file has wrong format\n");
    }
  }
  
  fin.close();
  nSnp_ = chrom_.size();
  tranLogDiff_.resize(nSnp_);
  tranLogSame_.resize(nSnp_);
  tranLogPlus_.resize(nSnp_);
}

int Parameter::validateLine(const string& line) const {

  istringstream sin(line);
  int count = 0;
  string token;
  while(sin >> token) {
    ++ count;
  }
  return count;
}

void Parameter::renewTran(uint32_t nHap) {

  nHap_ = nHap;
  for (size_t i = 1; i < nSnp_; i++) {
    double rho;
    if (gmap_[i] - gmap_[i - 1] < 1e-12) {
      rho = 4 * tranNe_ * 1e-12;
    } else {
      rho = 4 * tranNe_ * (gmap_[i] - gmap_[i - 1]);
    }

    double eRhoN = exp(-rho / nHap_);
    tranLogDiff_[i] = smoothMu_ * log((1 - eRhoN) / nHap_);
    tranLogSame_[i] = smoothMu_ * log(eRhoN + (1 - eRhoN) / nHap_);
    tranLogPlus_[i] = log(1 - exp(tranLogDiff_[i] - tranLogSame_[i])) + 
                     tranLogSame_[i];
  }
}
