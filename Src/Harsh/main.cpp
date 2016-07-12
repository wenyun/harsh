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

