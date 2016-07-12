#pragma once

#include <ctime>
#include <string>
#include <utility>
#include <vector>

using std::string;
using std::vector;
using std::pair;

class Read;
class Haplotype;
class Parameter;

class ReadIterator {

public:
  ReadIterator(const Read& read, const uint8_t readEnd, const uint8_t offset);
  pair<uint32_t, bool> operator*() const;
  void operator++();
  bool operator!=(const ReadIterator& iter) const;
private:
  const Read& read_;
  uint8_t readEnd_;
  uint8_t offset_;
};

class Read {
  
  friend class ReadIterator;

public:
  Read(const string& readLine, const Parameter& param);
  uint32_t countMatch(const Haplotype& hap) const;
  bool sampleFlip(const Haplotype& hap, const Parameter& param);
  bool sampleFlip(const Haplotype& hap1, 
                  const Haplotype& hap2,
                  const Parameter& param);
  
  typedef ReadIterator Iterator;
  Iterator begin() const;
  Iterator end() const;

  bool getFlip() const ;
  uint32_t size() const;
  void print() const;
  static clock_t getTime();
  static uint64_t getCount();
private:
  int validateReadEnd(string& readEnd) const;

  bool flip_;
  uint32_t size_;
  vector<uint32_t> position_;
  vector<uint16_t> startIndex_;
  vector<bool> code_;

  static clock_t time_;
  static uint64_t count_;
};
