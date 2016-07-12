#include "Reference.h"

#include "Haplotype.h"
#include "HarshIO.h"
#include "HarshUtil.h"
#include "Parameter.h"

#include <cfloat>
#include <cmath>
#include <ctime>
#include <ctime>
#include <fstream>
#include <string>

using std::ifstream;
using std::string;

extern char line[1000];

RefPath::RefPath(const Reference& reference) :
  reference_(reference),
  path_(reference.getSnpNo()) {}

bool RefPath::operator[](const uint32_t index) const {
  return reference_.getAllele(path_[index], index);
}

uint16_t& RefPath::operator()(const uint32_t index) {
  return path_[index];
}


int RefPath::size() const {
  return path_.size();
}

void RefPath::print() const {
  printf("Path\t");
  for (int i = 0; i < reference_.getSnpNo(); i++) {
    printf("%d_%d\t", path_[i], (*this)[i]);
  }
  printf("\n");
}

Reference::Reference()
  : nHap_(0),
    nSnp_(0),
    mode_(HapFirst),
    hap_(0) {}

void Reference::initFromFile(const char* refFile, const Parameter& param) {

  ifstream fin(refFile);
  string refLine;
  int lc = 0;
  while(getline(fin, refLine)) {
    int nBit = validateLine(refLine);
    if (nBit > 0) {
      ++ nSnp_;
      if (nHap_ == 0) {
        nHap_ = nBit;
      }
      if (nHap_ != nBit) {
        sprintf(line, "Reference file wrong format\n");
        harshErrorExit(line);
      }
      for (size_t i = 0; i < refLine.size(); i++) {
        if (refLine[i] == param.getMinorAllele(lc)) {
          hap_.push_back(true);
        } else if (refLine[i] == param.getMajorAllele(lc)) {
          hap_.push_back(false);
        } else {
          sprintf(line, 
                  "Reference file %s constains unexpected character %c for %d-th SNP (%c or %c expected)\n",
                  refFile, 
                  refLine[i],
                  lc + 1,
                  param.getMinorAllele(i),
                  param.getMajorAllele(i)
                  );
          harshErrorExit(line);
        }
      }
    }
    ++lc;
  }
  fin.close();
}

void Reference::samplePath(RefPath& path,
                           const Haplotype& haplotype,
                           const Parameter& param) const {
  
  clock_t begin = clock();

  vector<double> v(getHapNo() * getSnpNo());
  vector<double> vsum(getSnpNo());
  forward(v, vsum, haplotype, param);
  backwardSampling(path, v, vsum, param);

  clock_t end = clock();
  time_ += end - begin;
  ++ count_;
}

void plusTransition(vector<double>& tmp,
                    const double extra,
                    const int j,
                    const int k,
                    const int hapNo) {
  int i;
  for (i = 0; i < hapNo; i++) {
    tmp[k + i * hapNo] += extra;
    tmp[i + j * hapNo] += extra;
  }
}

int findMax(const vector<double>& score, 
            const int i, 
            const int hapNo, 
            double* maxScore) {
  int jk;

  double tmpScore = -DBL_MAX;
  int tmpIndex;  
  for (jk = 0; jk < hapNo * hapNo; jk++) {
    double scorejk = score[jk + i * hapNo * hapNo];
    if (scorejk > tmpScore) {
      tmpScore = scorejk;
      tmpIndex = jk;
    }
  }
  
  maxScore[0] = tmpScore;
  return tmpIndex;
}

void Reference::maxPath(RefPath& maxPath1,
                        RefPath& maxPath2,
                        const vector<pair<uint16_t, uint16_t> >& count,
                        const Parameter& param) const {

  Parameter maxParam(param);
  maxParam.renewTran(maxParam.getSampleSize());

  int hapNo = maxParam.getSampleSize();
  int snpNo = getSnpNo();
  int hapNo2 = hapNo * hapNo;
  double maxScore;
  vector<double> v(hapNo2 * 2, 0.0);
  vector<double> tmp(hapNo2, 0.0);
  vector<uint16_t> u(hapNo2 * snpNo, 0);

  vector<int> hapIdx;
  nSamplek(hapIdx, getHapNo(), maxParam.getSampleSize());

  plusEmission(v, 0, count[0], hapIdx, maxParam);
  
  for (int i = 1; i < snpNo; i++) {
    
    if (i % 1000 == 0) {
      printf(".");
      fflush(stdout);
    }

    for (int j = 0; j < hapNo; j++) {
      for (int k = 0; k < hapNo; k++) {
        tmp[k + j * hapNo] = v[k + j * hapNo + ((i - 1) % 2) * hapNo2] +
                             2 * maxParam.getTranLogDiff(i);
      }
    }
    
    for (int j = 0; j < hapNo; j++) {
      for (int k = 0; k < hapNo; k++) {
        plusTransition(tmp, 
                       maxParam.getTranLogSame(i) -
                       maxParam.getTranLogDiff(i),
                       j, 
                       k,
                       hapNo);
       
        u[k + j * hapNo + i * hapNo2] = findMax(tmp, 0, hapNo, &maxScore);
        v[k + j * hapNo + (i % 2) * hapNo2] = maxScore;

        plusTransition(tmp, 
                       maxParam.getTranLogDiff(i) -
                       maxParam.getTranLogSame(i),
                       j, 
                       k,
                       hapNo);
      }
    } 
    plusEmission(v, i, count[i], hapIdx, maxParam);
  }
  printf("\n");

  int maxJK = findMax(v, (snpNo - 1) % 2, hapNo, &maxScore);
  int maxK = maxJK % hapNo;
  int maxJ = maxJK / hapNo;
  maxPath1(snpNo - 1) = hapIdx[maxJ];
  maxPath2(snpNo - 1) = hapIdx[maxK];
  
  for (int i = snpNo - 1; i > 0; i--) {
    maxJK = u[maxK + maxJ * hapNo + i * hapNo2];
    maxK = maxJK % hapNo;
    maxJ = maxJK / hapNo;
    maxPath1(i - 1) = hapIdx[maxJ];
    maxPath2(i - 1) = hapIdx[maxK];
  }
}

void Reference::plusEmission(vector<double>& v,
                             const int i,
                             const pair<uint16_t, uint16_t>& count,
                             const vector<int>& hapIdx,
                             const Parameter& param) const {

  int hapNo = param.getSampleSize();
  int hapNo2 = hapNo * hapNo;

  int total = count.first + count.second;
  double prob0 = log(nChoosek(total, count.first)) +
                 count.first * param.getEmitLogRight() +
                 count.second * param.getEmitLogError();
  double prob1 = log(nChoosek(total, count.first)) +
                 total * log(0.5);
  double prob2 = log(nChoosek(total, count.first)) +
                 count.first * param.getEmitLogError() +
                 count.second * param.getEmitLogRight();
  
  for (int j = 0; j < hapNo; j++) {
    bool alleleJ = getAllele(hapIdx[j], i);
    for (int k = 0; k < hapNo; k++) {
      bool alleleK = getAllele(hapIdx[k], i);
      switch((int) (alleleJ + alleleK)) {
        case 0:
          v[k + j * hapNo + (i % 2) * hapNo2] += prob0;
          break;
        case 1:
          v[k + j * hapNo + (i % 2) * hapNo2] += prob1;
          break;
        case 2:
          v[k + j * hapNo + (i % 2) * hapNo2] += prob2;
          break;
        default:
          break;
      }
    }
  }
}

bool Reference::getAllele(const uint16_t hap, const uint32_t snp) const {
  if (mode_ == SnpFirst) {
    return hap_[snp * nHap_ + hap];
  } else if (mode_ == HapFirst) {
    return hap_[hap * nSnp_ + snp];
  }
  harshErrorExit("Bug found in Reference::getAllele");
  return true;
}

int Reference::getHapNo() const {
  return nHap_;
}

int Reference::getSnpNo() const {
  return nSnp_;
}

double Reference::forwardProb(const Haplotype& haplotype,
                              const Parameter& param) const {
   
  vector<double> v(getHapNo() * getSnpNo());
  vector<double> vsum(getSnpNo());
  forward(v, vsum, haplotype, param);
  return vsum[getSnpNo() - 1];
}

void Reference::forward(vector<double>& v, 
                        vector<double>& vsum,
                        const Haplotype& haplotype,
                        const Parameter& param) const {
  
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = 0.0;
  }
  
  for (size_t i = 0; i < vsum.size(); i++) {
    vsum[i] = 0.0;
  }

  for (int j = 0; j < getHapNo(); j++) {
    if (getAllele(j, 0) == haplotype[0]) {
      v[j] = param.getEmitLogRight();
    } else {
      v[j] = param.getEmitLogError();
    }
  }

  vsum[0] = logsumArray(v, 0, getHapNo());
  
  for (int i = 1; i < getSnpNo(); i++) {
    double tranBase = param.getTranLogDiff(i) + vsum[i - 1];
    for (int j = 0; j < getHapNo(); j++) {
      v[j + i * getHapNo()] = 
        logsum(v[j + (i - 1) * getHapNo()] + param.getTranLogPlus(i), 
               tranBase);
    }
    for (int j = 0; j < getHapNo(); j++) {
      if (getAllele(j, i) == haplotype[i]) {
        v[j + i * getHapNo()] += param.getEmitLogRight();
      } else {
        v[j + i * getHapNo()] += param.getEmitLogError();
      }
    }
  
    vsum[i] = logsumArray(v, i * getHapNo(), getHapNo());
  }
}

void Reference::backwardSampling(RefPath& path,
                                 const vector<double>& v,
                                 const vector<double>& vsum,
                                 const Parameter& param) const {

  vector<double> prob(getHapNo());
  vector<double> tran(getHapNo());

  for (int k = 0; k < getHapNo(); k++) {
    prob[k] = exp(v[k + (getSnpNo() - 1) * getHapNo()] - vsum[getSnpNo() - 1]);
  }
  path(getSnpNo() - 1) = mnrnd(prob) - 1;
  
  for (int i = getSnpNo() - 2; i >= 0; i--) {
    for (int j = 0; j < getHapNo(); j++) {
      tran[j] = v[j + i * getHapNo()] + param.getTranLogDiff(i + 1);
    }

    tran[path(i + 1)] = v[path(i + 1) + i * getHapNo()] + 
                        param.getTranLogSame(i + 1);

    double tranSum = 
      logsum(vsum[i] + param.getTranLogDiff(i + 1),
             v[path(i + 1) + i * getHapNo()] + param.getTranLogPlus(i + 1));
    
    for (int k = 0; k < getHapNo(); k++) {
      prob[k] = exp(tran[k] - tranSum);
    }
    path(i) = mnrnd(prob) - 1;
  }
}

int Reference::validateLine(string& line) const {
  
  size_t nBit = 0;
  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] != '-' && 
        line[i] != ' ' && 
        line[i] != '\t' && 
        line[i] != '\r' && 
        line[i] != '\n') {
      ++ nBit;
    }
  }

  if (nBit != line.size()) {
    string clean(nBit, ' ');
    int j = 0;
    for (size_t i = 0; i < line.size(); i++) {
      if (line[i] != '-' && 
          line[i] != ' ' && 
          line[i] != '\t' && 
          line[i] != '\r' && 
          line[i] != '\n') {
        clean[j++] = line[i];
      }
    }
    line = clean;
  }

  return nBit;
}

clock_t Reference::getTime() {
  return time_;
}

uint64_t Reference::getCount() {
  return count_;
}

clock_t Reference::time_ = 0;
uint64_t Reference::count_ = 0;
