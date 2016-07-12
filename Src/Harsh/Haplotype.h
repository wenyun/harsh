#pragma once

#include <ctime>
#include <string>
#include <vector>
#include <utility>

using std::string;
using std::vector;
using std::pair;

class Read;
class RefPath;
class Parameter;

class Haplotype {

public:
  Haplotype(const int length);
  Haplotype(const vector<bool>& haplotype);

  void copyInit(const RefPath& path);
  enum Mode{NonFlipOnly, FlipOnly, AllWithFlipRev};  
  bool assignRead(const Read& read, const Mode mode);
  bool removeRead(const Read& read, const Mode mode);
  bool flipRead(const Read& read, const Mode mode);
  void sample(const RefPath& path, const Parameter& param);
  void sample(const Parameter& param);

  int getLength() const;
  int getCount(const uint32_t index, const bool allele) const;

  bool operator[](const uint32_t index) const;
  void print() const;

  static clock_t getTime();
  static uint64_t getCount();

private:
  void randomInit();
  void appendRead(const Read& read, const bool flip);
  void popRead(const Read& read, const bool flip);

  uint32_t length_;
  vector<pair<uint16_t, uint16_t> > readCount_;
  vector<bool> haplotype_;

  static clock_t time_;
  static uint64_t count_;
};
