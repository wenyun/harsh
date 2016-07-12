#pragma once

#include <ctime>
#include <vector>
#include <utility>
#include <string>

using std::vector;
using std::pair;
using std::string;

class Haplotype;
class Parameter;
class Reference;

class RefPath {

public:
  RefPath(const Reference& reference);
  bool operator[](const uint32_t index) const;
  uint16_t& operator()(const uint32_t index);
  int size() const;
  void print() const;
private:
  const Reference& reference_;
  vector<uint16_t> path_;  
};

class Reference {
  
  friend class RefPath;

public:
  
  Reference();
  void initFromFile(const char* refFile, const Parameter& param);

  typedef RefPath Path;
  void samplePath(RefPath& path,
                  const Haplotype& hap,
                  const Parameter& param) const;
  
  void maxPath(RefPath& path1,
               RefPath& path2,
               const vector<pair<uint16_t, uint16_t> >& count,
               const Parameter& param) const;

  double forwardProb(const Haplotype& haplotype, 
                     const Parameter& param) const;

  int getHapNo() const;
  int getSnpNo() const;
  
  static clock_t getTime();
  static uint64_t getCount();

private:

  bool getAllele(const uint16_t hap, const uint32_t snp) const;
  
  void forward(vector<double>& v, 
               vector<double>& vsum,
               const Haplotype& haplotype,
               const Parameter& param) const;

  void backwardSampling(RefPath& path,
                        const vector<double>& v,
                        const vector<double>& vsum,
                        const Haplotype& haplotype,
                        const Parameter& param) const;

  void plusEmission(vector<double>& v,
                    const int i,
                    const pair<uint16_t, uint16_t>& count,
                    const vector<int>& hapIdx,
                    const Parameter& param) const;

  int validateLine(string& line) const;

  enum ModeType{SnpFirst, HapFirst};
  int nHap_;
  int nSnp_;
  ModeType mode_;
  vector<bool> hap_;

  static clock_t time_;
  static uint64_t count_;
};
