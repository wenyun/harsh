#include "Harsh.h"

#include "Haplotype.h"
#include "HarshIO.h"
#include "Parameter.h"
#include "Read.h"
#include "Reference.h"

#include <algorithm> 
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <climits>
#include <cfloat>

using std::make_pair;
using std::min;
using std::max;

extern char line[1000];

void exitWithHelp() {
  printf(
  "Usage: harsh [options]\n"
  "options: \n"
  "--sfile seq_file : sequencing read file\n"
  "--rfile ref_file : reference haplotype file\n"
  "--mfile map_file : SNP map file\n"
  "--output out_file : output haplotype file\n"
  "-n iteration: number of sampling iterations (default 1e3)\n"
  "-u mu : smoothing parameter for sampling (default 1)\n"
  "-e epsilon : sequencing error rate (default 0.01)\n"
  "-w omega : mismatch error rate between reference and donar haplotype (default 2e-3)\n"
  "-r rho : population constant for computing transitions (default 1e-3)\n"
  "-v verbose : verbose level (default 1)\n"
  );
  exit(1);
}

void printVersionInformation() {
  printf(
  "@----------------------------------------------------------@\n"
  "|         HARSH!       |      v0.21      |    8/Nov/201 3  |\n"
  "|----------------------------------------------------------|\n"
  "|                    (C) 2013 Wen-Yun Yang                 |\n"
  "|----------------------------------------------------------|\n"
  "|  For documentation, citation & bug-report instructions:  |\n"
  "|            http://genetics.cs.ucla.edu/harsh             |\n"
  "@----------------------------------------------------------@\n"
  "\n"
  );
}

int main(int argc, char** argv) {

  printVersionInformation();
  srand(time(NULL));

  Parameter param;
  parseInputParameter(argc, argv, param);
  param.initFromMapFile();

  vector<Read> read;
  readSeqFile(param.getSeqFile(), read, param);

  Reference reference;
  if (param.getHasRef()) {
    reference.initFromFile(param.getRefFile(), param);
    param.setSampleSize(min(reference.getHapNo(), param.getSampleSize()));
    param.renewTran(reference.getHapNo());
  }
  
  Haplotype haplotype1(param.getNSnp());
  Haplotype haplotype2(param.getNSnp());
  
  if (param.getHasRef()) {
    /* haplotype inference */
    printInputInformation(read, reference, param);
    harshMessage("Haplotype inference with reference start...", 1, param);
    inferenceGibbs(haplotype1, haplotype2, read, reference, param);
    outputInference(haplotype1, haplotype2, param);
  } else {
    /* haplotype assembly */
    printInputInformation(read, param);
    harshMessage("Haplotype assembly start...", 1, param);
    assemblyGibbs(haplotype1, read, param);
    printf("start output\n");
    outputAssembly(haplotype1, param);
  }

  printf("sampling time (read %llu %f, haplotype %llu %f, reference %llu %f)\n",
         Read::getCount(),
         ((double) Read::getTime()) / CLOCKS_PER_SEC,
         Haplotype::getCount(),
         ((double) Haplotype::getTime()) / CLOCKS_PER_SEC,
         Reference::getCount(),
         ((double) Reference::getTime()) / CLOCKS_PER_SEC);
}

void printInputInformation(const vector<Read>& read,
                           const Reference& reference,
                           const Parameter& param) {
  
  sprintf(line,
  "Options in effect:\n"
  "    --sfile %s\n"
  "    --rfile %s\n"
  "    --mfile %s\n"
  "    --output %s\n"
  "    -n %d\n"
  "    -u %f\n"
  "    -e %f\n"
  "    -w %f\n"
  "    -r %f\n"
  "    -v %d\n",
  param.getSeqFile(),
  param.getRefFile(),
  param.getMapFile(),
  param.getOutFile(),
  param.getMaxIter(),
  param.getSmoothMu(),
  param.getSeqError(),
  param.getEmitError(),
  param.getTranNe(),
  param.getVerbose()
  );
  harshMessage(line, 0, param);

  int64_t nBit = 0;
  for (int i = 0; i < read.size(); i++) {
    nBit += read[i].size();
  }

  sprintf(line, "Number of read : [%d]", (int)read.size());
  harshMessage(line, 0, param);  
  sprintf(line, "Number of read bit : [%lld]", nBit);
  harshMessage(line, 0, param);  
  sprintf(line, "Number of SNP per read on average : [%.2f]", ((double) nBit) / read.size());
  harshMessage(line, 0, param);  
  sprintf(line, "Number of SNP : [%d]", param.getNSnp());
  harshMessage(line, 0, param);  
  sprintf(line, "Number of reference : [%d]", reference.getHapNo());
  harshMessage(line, 0, param);  
  harshMessage("", 0, param);  
}


void printInputInformation(const vector<Read>& read,
                           const Parameter& param) {
  
  sprintf(line,
  "Options in effect:\n"
  "    --sfile %s\n"
  "    --mfile %s\n"
  "    --output %s\n"
  "    -n %d\n"
  "    -u %f\n"
  "    -e %f\n"
  "    -w %f\n"
  "    -r %f\n"
  "    -v %d\n",
  param.getSeqFile(),
  param.getMapFile(),
  param.getOutFile(),
  param.getMaxIter(),
  param.getSmoothMu(),
  param.getSeqError(),
  param.getEmitError(),
  param.getTranNe(),
  param.getVerbose()
  );
  harshMessage(line, 0, param);

  int64_t nBit = 0;
  for (int i = 0; i < read.size(); i++) {
    nBit += read[i].size();
  }

  sprintf(line, "Number of read : [%d]", (int)read.size());
  harshMessage(line, 0, param);  
  sprintf(line, "Number of read bit : [%lld]", nBit);
  harshMessage(line, 0, param);  
  sprintf(line, "Number of SNP per read on average : [%.2f]", ((double) nBit) / read.size());
  harshMessage(line, 0, param);  
  sprintf(line, "Number of SNP : [%d]", param.getNSnp());
  harshMessage(line, 0, param);  
  harshMessage("", 0, param);  
}


void parseInputParameter(int argc,
                         char** argv,
                         Parameter& param) {
  int i;
  
  if (argc == 1) {
    exitWithHelp();
  }

  for(i = 1; i < argc; i++) {
    if(argv[i][0] != '-') {
      exitWithHelp();
    }
    if(++i >= argc) {
      exitWithHelp();
    }

    switch(argv[i-1][1]) {
      case '-':
        if(!strcmp(argv[i-1], "--sfile")) {
          param.setSeqFile(argv[i]);
        } else if(!strcmp(argv[i-1], "--rfile")) {
          param.setRefFile(argv[i]);
        } else if(!strcmp(argv[i-1], "--mfile")) {
          param.setMapFile(argv[i]);
        } else if(!strcmp(argv[i-1], "--output")) {
          param.setOutFile(argv[i]);
        } else {
          sprintf(line, "Unknown options: %s\n", argv[i-1]);
          harshWarning(line);
          exitWithHelp();
        }
        break;
      case 'n':
        param.setMaxIter(atoi(argv[i]));
        break;
      case 'm':
        param.setSubIter(atoi(argv[i]));
        break;        
      case 'u':
        param.setSmoothMu(atof(argv[i]));
        break;
      case 'e':
        param.setSeqError(atof(argv[i]));
        break;
      case 'w':
        param.setEmitError(atof(argv[i]));
        break;
      case 'r':
        param.setTranNe(atof(argv[i]));
        break;
      case 'v':
        param.setVerbose(atoi(argv[i]));
        break;
      case 's':
        param.setSampleSize(atoi(argv[i]));
        break;
      default:
        sprintf(line, "Unknown options: -%s\n", argv[i-1]);
        harshWarning(line);
        exitWithHelp();
    }
  }
}

void assemblyGibbs(Haplotype& haplotype,
                   vector<Read>& read,
                   const Parameter& param) {

  /* initialize read count */
  for (int i = 0; i < read.size(); i++) {
    haplotype.assignRead(read[i], Haplotype::AllWithFlipRev);
  }
 
  int minMismatch = INT_MAX;
  Haplotype optimal = haplotype;
  /* sampling iteration */
  for (int i = 0; i < param.getMaxIter(); i++) {
    haplotype.sample(param);
    for (int j = 0; j < read.size(); j++) {
      bool isChanged = read[j].sampleFlip(haplotype, param);
      if (isChanged) {
        haplotype.flipRead(read[j], Haplotype::AllWithFlipRev);
      }
    }
    
    int mismatch = assemblyMismatch(haplotype);
    sprintf(line, "no. of mismatch = %d", mismatch);
    harshMessage(line, 2, param);

    if (mismatch < minMismatch) {
      minMismatch = mismatch;
      optimal = haplotype;
    }
  }
  
  haplotype = optimal;

  sprintf(line, "Optimal number of mismatch : [%d]", minMismatch);
  harshMessage(line, 0, param);
}

int assemblyMismatch(Haplotype& haplotype) {
  int mismatch = 0;
  for (int i = 0; i < haplotype.getLength(); i++) {
    if (haplotype[i]) {
      mismatch += haplotype.getCount(i, 0);
    } else {
      mismatch += haplotype.getCount(i, 1);
    }
  }
  return mismatch;
}

void inferenceGibbs(Haplotype& haplotype1,
                    Haplotype& haplotype2,
                    vector<Read>& read,
                    const Reference& reference,
                    const Parameter& param) {

  if (param.getMaxIter() > 0) {
    harshMessage("counting read likelihood...", 2, param);

    /* make the read count likelihood */
    vector<pair<uint16_t, uint16_t> > count(reference.getSnpNo(), 
                                            make_pair<uint16_t, uint16_t>(0, 0));
    for (int i = 0; i < read.size(); i++) {
      for (Read::Iterator iter = read[i].begin(); 
           iter != read[i].end(); 
           ++iter) {
        uint32_t index = (*iter).first;
        bool bit = (*iter).second;
        if (!bit) {
          ++ count[index].first;
        } else {
          ++ count[index].second;
        }
      }
    }
    
    harshMessage("warm start...", 1, param);
    
    /* find the best haplotype imputation */
    RefPath maxPath1(reference);
    RefPath maxPath2(reference);
    reference.maxPath(maxPath1, maxPath2, count, param);
    haplotype1.copyInit(maxPath1);
    haplotype2.copyInit(maxPath2);
    
    /* sample read assignment */
    for (int i = 0; i < read.size(); i++) {
      read[i].sampleFlip(haplotype1, haplotype2, param);
    }
  }

  double likelihood = 
    inferenceLikelihood(haplotype1, haplotype2, read, reference, param);
  sprintf(line, "Start likelihood = %f\n", likelihood);
  harshMessage(line, 2, param);  

  /* initialize haplotype read counts */
  for (int i = 0; i < read.size(); i++) {
    if (!read[i].getFlip()) {
      haplotype1.assignRead(read[i], Haplotype::NonFlipOnly);
    } else {
      haplotype2.assignRead(read[i], Haplotype::FlipOnly);
    }
  }

  double maxLikelihood = -DBL_MAX;
  Haplotype optimal1 = haplotype1;
  Haplotype optimal2 = haplotype2;

  RefPath path1(reference);
  RefPath path2(reference);
  /* sampling iteration */
  for (int i = 0; i < abs(param.getMaxIter()); i++) {
    reference.samplePath(path1, haplotype1, param);
    reference.samplePath(path2, haplotype2, param);

    for (int k = 0; k < param.getSubIter(); k++) {
      haplotype1.sample(path1, param);
      haplotype2.sample(path2, param);
      for (int j = 0; j < read.size(); j++) {
        bool isChanged = read[j].sampleFlip(haplotype1, haplotype2, param);
        if (isChanged) {
          if (read[j].getFlip()) {
            haplotype1.removeRead(read[j], Haplotype::NonFlipOnly);
            haplotype2.assignRead(read[j], Haplotype::FlipOnly);
          } else {
            haplotype1.assignRead(read[j], Haplotype::NonFlipOnly);
            haplotype2.removeRead(read[j], Haplotype::FlipOnly);
          }
        }
      }
    }
    
    double likelihood = 
      inferenceLikelihood(haplotype1, haplotype2, read, reference, param);
    sprintf(line, "likelihood = %f (%.1f%% done)", 
            likelihood, 
            ((double) i) / abs(param.getMaxIter())* 100);
    harshMessage(line, 1, param);

    if (likelihood > maxLikelihood) {
      maxLikelihood = likelihood;
      optimal1 = haplotype1;
      optimal2 = haplotype2;
    }
  }

  sprintf(line, "max likelihood = [%.5f]", maxLikelihood);
  harshMessage(line, 1, param);
  
  haplotype1 = optimal1;
  haplotype2 = optimal2;
}

double inferenceLikelihood(Haplotype& haplotype1,
                           Haplotype& haplotype2,
                           vector<Read>& read,
                           const Reference& reference,
                           const Parameter& param) {

  double logLikelihoodRead = 0.0;
  double logLikelihoodRef = 0.0;
  uint32_t mismatchEither = 0;
  uint32_t mismatchBoth = 0;
  uint32_t matchEither = 0;
  uint32_t matchBoth = 0;
  for (int i = 0; i < read.size(); i++) {
    uint32_t match1 = read[i].countMatch(haplotype1);
    uint32_t match2 = read[i].countMatch(haplotype2);
    uint32_t mismatch1 = read[i].size() - match1;
    uint32_t mismatch2 = read[i].size() - match2;

    mismatchEither += min(mismatch1, mismatch2);
    mismatchBoth += mismatch1 + mismatch2;
    matchEither += max(match1, match2);
    matchBoth += match1 + match2;
  }

  logLikelihoodRead += param.getSeqLogRight() * matchBoth + 
                       param.getSeqLogError() * mismatchBoth;

  logLikelihoodRef += reference.forwardProb(haplotype1, param);
  logLikelihoodRef += reference.forwardProb(haplotype2, param);

  sprintf(line, "mismatch = %.5f [%d %d], reference = %.5f",
          logLikelihoodRead,
          mismatchBoth, 
          mismatchEither,
          logLikelihoodRef);
  harshMessage(line, 3, param);

  return logLikelihoodRead + logLikelihoodRef;
}
