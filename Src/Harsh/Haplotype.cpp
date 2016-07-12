#include "Haplotype.h"

#include "HarshIO.h"
#include "Parameter.h"
#include "Read.h"
#include "Reference.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

Haplotype::Haplotype(const int length)
  : length_(length),
    readCount_(length),
    haplotype_(length) {
  
  randomInit();
}

Haplotype::Haplotype(const vector<bool>& haplotype)
  : length_(haplotype.size()),
    readCount_(haplotype.size()),
    haplotype_(haplotype) {}

void Haplotype::copyInit(const RefPath& path) {
  for (int i = 0; i < length_; i++) {
    haplotype_[i] = path[i];
  }
}

bool Haplotype::assignRead(const Read& read, const Mode mode) {
  
  bool done = false;
  switch (mode) {
    case NonFlipOnly:
      if (!read.getFlip()) {
        appendRead(read, false);
        done = true;
      } 
      break;
    case FlipOnly:
      if (read.getFlip()) {
        appendRead(read, false);
        done = true;
      } 
      break;
    case AllWithFlipRev:
      appendRead(read, read.getFlip());
      done = true;
      break;
    default:
      harshErrorExit("You should not get to here\n");
  }
  return done;
}

bool Haplotype::removeRead(const Read& read, const Mode mode) {
  
  bool done = false;
  switch (mode) {
    case NonFlipOnly:
      if (read.getFlip()) {
        popRead(read, false);
        done = true;
      }  
      break;
    case FlipOnly:
      if (!read.getFlip()) {
        popRead(read, false);
        done = true;
      } 
      break;
    default:
      harshErrorExit("You should not get here\n");
  }
  return done;
}

bool Haplotype::flipRead(const Read& read, const Mode mode) {

  bool done = false;
  if (mode == AllWithFlipRev) {
    if (read.getFlip()) {
      popRead(read, false);
      appendRead(read, true);
    } else {
      popRead(read, true);
      appendRead(read, false);
    }
    done = true;
  }
  return done;
}

void Haplotype::sample(const RefPath& path, const Parameter& param) {
  
  clock_t begin = clock();

  for (int i = 0; i < getLength(); i++) {
    uint16_t count0 = readCount_[i].first;
    uint16_t count1 = readCount_[i].second;
    double weight0 = exp(param.getSmoothMu() * 
                         (param.getSeqLogRight() * count0 +
                          param.getSeqLogError() * count1 +
                          param.getEmitLogRight() * (path[i] == 0) +
                          param.getEmitLogError() * (path[i] == 1))
                        );
    double weight1 = exp(param.getSmoothMu() * 
                         (param.getSeqLogError() * count0 +
                          param.getSeqLogRight() * count1 +
                          param.getEmitLogError() * (path[i] == 0) +
                          param.getEmitLogRight() * (path[i] == 1))
                        );
    double prob0 = weight0 / (weight0 + weight1);


    if (((double) rand()) / RAND_MAX < prob0) {
      haplotype_[i] = false;
    } else {
      haplotype_[i] = true;
    }
  }

  clock_t end = clock();
  time_ += end - begin;
  ++ count_;
}

void Haplotype::sample(const Parameter& param) {

  clock_t begin = clock();
  
  for (int i = 0; i < getLength(); i++) {
    uint16_t count0 = readCount_[i].first;
    uint16_t count1 = readCount_[i].second;
    double weight0 = exp(param.getSmoothMu() * 
                         (param.getSeqLogRight() * count0 +
                          param.getSeqLogError() * count1)
                        );
    double weight1 = exp(param.getSmoothMu() * 
                         (param.getSeqLogError() * count0 +
                          param.getSeqLogRight() * count1)
                        );
    double prob0 = weight0 / (weight0 + weight1);

    if (((double) rand()) / RAND_MAX < prob0) {
      haplotype_[i] = false;
    } else {
      haplotype_[i] = true;
    }
  }

  clock_t end = clock();
  time_ += end - begin;  
  ++ count_;
}

int Haplotype::getLength() const {
  return length_;
}

int Haplotype::getCount(const uint32_t index, const bool allele) const {
  if (!allele) {
    return readCount_[index].first;
  } else {
    return readCount_[index].second;
  }
}


bool Haplotype::operator[](const uint32_t index) const {
  return haplotype_[index];
}

void Haplotype::randomInit() {
  for (int i = 0; i < length_; i++) {
    if (rand() % 2 == 0) {
      haplotype_[i] = false;
    } else {
      haplotype_[i] = true;
    }
  }
}

void Haplotype::appendRead(const Read& read, const bool flip) {
  for (Read::Iterator iter = read.begin(); iter != read.end(); ++iter) {
    uint32_t index = (*iter).first;
    bool bit = (*iter).second;
    if (bit ^ flip) {
      ++ (readCount_[index].second);
    } else {
      ++ (readCount_[index].first);
    }
  }
}

void Haplotype::popRead(const Read& read, const bool flip) {
  for (Read::Iterator iter = read.begin(); iter != read.end(); ++iter) {
    uint32_t index = (*iter).first;
    bool bit = (*iter).second;
    if (bit ^ flip) {
      -- (readCount_[index].second);
    } else {
      -- (readCount_[index].first);
    }
  }
}

void Haplotype::print() const {
  printf("hap\t");
  for (int i = 0; i < getLength(); i++) {
    printf("%d\t", (*this)[i]);
  }
  printf("\n");
  printf("cnt0\t");
  for (int i = 0; i < getLength(); i++) {
    printf("%d\t", (*this).getCount(i, 0));
  }
  printf("\n");
  printf("cnt1\t");
  for (int i = 0; i < getLength(); i++) {
    printf("%d\t", (*this).getCount(i, 1));
  }
  printf("\n");
}

clock_t Haplotype::getTime() {
  return time_;
}

uint64_t Haplotype::getCount() {
  return count_;
}

clock_t Haplotype::time_ = 0;
uint64_t Haplotype::count_ = 0;
