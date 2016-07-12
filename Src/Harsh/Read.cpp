#include "Read.h"

#include "Haplotype.h"
#include "Parameter.h"
#include "HarshIO.h"

#include <cmath>
#include <ctime>
#include <sstream>

using std::istringstream;

extern char line[1000];

ReadIterator::ReadIterator(const Read& read,
                           uint8_t readEnd,
                           uint8_t offset) : 
  read_(read),
  readEnd_(readEnd),
  offset_(offset) {}

pair<uint32_t, bool> ReadIterator::operator*() const {  
  return std::make_pair<uint32_t, bool>(
      read_.position_[readEnd_] + offset_,
      read_.code_[read_.startIndex_[readEnd_] + offset_ * 2 + 1]);
}

void ReadIterator::operator++() {
  do {
    if (read_.startIndex_[readEnd_] + offset_ * 2 + 2 == read_.code_.size()) {
      ++ readEnd_;
      offset_ = 0;
      break;
    } else if (read_.startIndex_[readEnd_] + offset_ * 2 + 2 == 
               read_.startIndex_[readEnd_ + 1]) {
      ++ readEnd_;
      offset_ = 0;
    } else {
      ++ offset_ ;
    }
  } while (read_.code_[read_.startIndex_[readEnd_] + offset_ * 2] == 0); 

}

bool ReadIterator::operator!=(const ReadIterator& iter) const {
  return (iter.readEnd_ != readEnd_) || (iter.offset_ != offset_);
}
  
Read::Read(const string& readLine, const Parameter& param)
  : flip_(0),
    size_(0),
    position_(0),
    startIndex_(0),
    code_(0) {
  
  flip_ = rand() % 2;
  uint32_t position = 0;
  string readEnd;
  uint16_t index = 0;

  istringstream iss(readLine);
  while (iss.good()) {
    iss >> position >> readEnd;
    /* change to index starting 0*/
    -- position;
  
    int beginDumpSize = validateReadEnd(readEnd);
    if (readEnd.size() == 0) {
      continue;
    } else {
      position += beginDumpSize;
    }
    
    position_.push_back(position);
    startIndex_.push_back(index);
    index += readEnd.size() * 2;

    for (size_t i = 0; i < readEnd.size(); i++) {
      if (readEnd[i] == '-') {
        code_.push_back(false);
        code_.push_back(false);
      } else if (readEnd[i] == param.getMinorAllele(position + i)) {
        code_.push_back(true);
        code_.push_back(true);
        ++ size_;
      } else if (readEnd[i] == param.getMajorAllele(position + i)) {
        code_.push_back(true);
        code_.push_back(false);
        ++ size_;
      } else {
        printf("%s %s\n", readLine.c_str(), readEnd.c_str());
        sprintf(line, 
                "[%s] - Inconsistent allele %c at position %ld, given %c, should be %c or %c",
                readLine.c_str(),
                readEnd[i],
                position + 1 + i,
                readEnd[i],
                param.getMinorAllele(position + i),
                param.getMajorAllele(position + i));
        harshWarning(line);
      }
    }
  }
}

uint32_t Read::countMatch(const Haplotype& hap) const {
  uint32_t match = 0;
  for (Iterator iter = begin(); iter != end(); ++iter) {
    uint32_t index = (*iter).first;
    bool bit = (*iter).second;
    if (hap[index] == bit) {
      ++ match;
    }
  }
  return match;
}

bool Read::sampleFlip(const Haplotype& hap, const Parameter& param) {

  clock_t begin = clock();

  uint32_t match = countMatch(hap);
  uint32_t mismatch = size() - match;
   
  double weight0 = exp(param.getSmoothMu() * 
                       (param.getSeqLogRight() * match +
                        param.getSeqLogError() * mismatch)
                      );
  double weight1 = exp(param.getSmoothMu() * 
                       (param.getSeqLogError() * match +
                        param.getSeqLogRight() * mismatch)
                      );

  double prob0 = weight0 / (weight0 + weight1);

  bool flipOld = flip_;
  if (((double) rand()) / RAND_MAX < prob0) {
    flip_ = false;
  } else {
    flip_ = true;
  }

  clock_t end = clock();;
  time_ += end - begin;
  ++ count_;

  if (flipOld == flip_) {
    return false;
  } else {
    return true;
  }
}

bool Read::sampleFlip(const Haplotype& hap1,
                      const Haplotype& hap2,
                      const Parameter& param) {

  clock_t begin = clock();
  
  uint32_t match1 = countMatch(hap1);
  uint32_t mismatch1 = size() - match1;
  uint32_t match2 = countMatch(hap2);
  uint32_t mismatch2 = size() - match2;

  double weight1 = exp(param.getSmoothMu() * 
                       (param.getSeqLogRight() * match1 +
                        param.getSeqLogError() * mismatch1)
                      );
  double weight2 = exp(param.getSmoothMu() * 
                       (param.getSeqLogRight() * match2 +
                        param.getSeqLogError() * mismatch2)
                      );

  double prob1 = weight1 / (weight1 + weight2);
  
  bool flipOld = flip_;
  if (((double) rand()) / RAND_MAX < prob1) {
    flip_ = false;
  } else {
    flip_ = true;
  }

  clock_t end = clock();
  time_ += end - begin;
  ++ count_;

  if (flipOld == flip_) {
    return false;
  } else {
    return true;
  }
}

Read::Iterator Read::begin() const {
  Iterator beginIter(*this, 0, 0);
  return beginIter;
}

Read::Iterator Read::end() const {
  Iterator endIter(*this, position_.size(), 0);
  return endIter;
}

bool Read::getFlip() const {
  return flip_;
}

uint32_t Read::size() const {
  return size_;
}

int Read::validateReadEnd(string& readEnd) const {
  /* remove beginning missing data */
  int beginDumpSize = 0;
  for (size_t i = 0; i < readEnd.size(); ++i) {
    if (readEnd[i] != '-') {
      break;
    } else {
      ++ beginDumpSize;
    } 
  }
  
  if (beginDumpSize > 0) {
    readEnd.erase(0, beginDumpSize);
  }

  /* remove tail missing data */
  int tailDumpSize = 0;
  for (int i = readEnd.size() - 1; i >= 0; --i) {
    if (readEnd[i] != '-') {
      break;
    } else {
      ++ tailDumpSize;
    } 
  }

  if (tailDumpSize > 0) {
    readEnd.erase(readEnd.size() - tailDumpSize, tailDumpSize);
  }
  
  return beginDumpSize;
}

void Read::print() const {
  
  int maxPos = 0;
  for (Iterator iter = begin(); iter != end(); ++iter) {
    int pos = (*iter).first;
    if (pos > maxPos) {
      maxPos = pos;
    }
  }

  vector<char> all(maxPos + 1, '-');
  for (Iterator iter = begin(); iter != end(); ++iter) {
    int pos = (*iter).first;
    bool bit = (*iter).second;
    if (bit) {
      all[pos] = '1';
    } else {
      all[pos] = '0';
    }
  }
    
  printf("Read\t");
  for (size_t i = 0; i < all.size(); i++) {
    if (getFlip() && all[i] != '-') {
      printf("(%c)\t", all[i]);
    } else {
      printf("%c\t", all[i]);
    }
  }
  printf("\n");
}

clock_t Read::getTime() {
  return time_;
}

uint64_t Read::getCount() {
  return count_;
}

clock_t Read::time_ = 0;
uint64_t Read::count_ = 0;
