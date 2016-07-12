#include "HarshIO.h"

#include "Haplotype.h"
#include "Parameter.h"
#include "Read.h"

#include <fstream>
#include <string>
#include <vector>

using std::ifstream;
using std::ofstream;
using std::endl;
using std::string;
using std::vector;

char line[1000];

void readSeqFile(const char* seqFile, 
                 vector<Read>& read,
                 const Parameter& param) {
  ifstream fin(seqFile);
  string seqLine;

  while(getline(fin, seqLine)) {
    if (!isEmpty(seqLine)) {
      Read readNew(seqLine, param);
      if (readNew.size() > 0) {
        read.push_back(readNew);
      }
    }
  }
  fin.close();
}

void outputAssembly(const Haplotype& haplotype,
                    const Parameter& param) {

  ofstream fout(param.getOutFile());
  for (int i = 0; i < haplotype.getLength(); i++) {
    fout << (haplotype[i] ? 
             param.getMinorAllele(i) : 
             param.getMajorAllele(i))
         << ((!haplotype[i]) ?
             param.getMinorAllele(i) : 
             param.getMajorAllele(i))
         << endl;
  }
  fout.close();
}

void outputInference(const Haplotype& haplotype1,
                     const Haplotype& haplotype2,
                     const Parameter& param) {

  ofstream fout(param.getOutFile());
  for (int i = 0; i < haplotype1.getLength(); i++) {
    fout << (haplotype1[i] ? 
             param.getMinorAllele(i) : 
             param.getMajorAllele(i))
         << (haplotype2[i] ? 
             param.getMinorAllele(i) : 
             param.getMajorAllele(i))
         << endl;
  }
  fout.close();
}

bool isEmpty(const string& line) {
  int nBit = 0;
  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] != '-' && line[i] != ' ' && line[i] != '\t') {
      ++ nBit;
    }
  }

  if (nBit == 0) {
    return true;
  } else {
    return false;
  }
}

void harshErrorExit(const char* msg) {
  fprintf(stderr, "FATAL => %s\n", msg);
  exit(1);
}

void harshWarning(const char* msg) {
  fprintf(stderr, "WARINING => %s\n", msg);
}

void harshMessage(const char* msg, 
                  const int level,
                  const Parameter& param) {
  if(level <= param.getVerbose()) {
    printf("%s\n", msg);
  }
}
