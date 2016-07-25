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
#include <iostream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
using std::make_pair;
using std::min;
using std::max;

extern char line[1000];

void printVersionInformation() {
    printf(
           "@----------------------------------------------------------@\n"
           "|         HARSH!       |      v0.3       |   25/Jul/2016   |\n"
           "|----------------------------------------------------------|\n"
           "|                    (C) 2016 Wen-Yun Yang                 |\n"
           "|----------------------------------------------------------|\n"
           "|  For documentation, citation & bug-report instructions:  |\n"
           "|            http://genetics.cs.ucla.edu/harsh             |\n"
           "@----------------------------------------------------------@\n"
           "\n"
           );
}

constexpr auto kSeqFile = "sfile";
constexpr auto kRefFile = "rfile";
constexpr auto kMapFile = "mfile";
constexpr auto kOutFile = "output";
constexpr auto kMaxIter = "n";
constexpr auto kSmoothMu = "u";
constexpr auto kSeqError = "e";
constexpr auto kEmitError = "w";
constexpr auto kTranNe = "r";
constexpr auto kVerbose = "v";

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
    for (size_t i = 0; i < read.size(); i++) {
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
    for (size_t i = 0; i < read.size(); i++) {
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

void parseInputParameter(const po::variables_map &vm,
                         Parameter& param) {
    
    if (vm.count(kSeqFile)) {
        param.setSeqFile(vm[kSeqFile].as<std::string>().c_str());
    }
    if (vm.count(kMapFile)) {
        param.setMapFile(vm[kMapFile].as<std::string>().c_str());
    }
    if (vm.count(kOutFile)) {
        param.setOutFile(vm[kOutFile].as<std::string>().c_str());
    }
    if (vm.count(kRefFile)) {
        param.setRefFile(vm[kRefFile].as<std::string>().c_str());
    }
    if (vm.count(kMaxIter)) {
        param.setMaxIter(vm[kMaxIter].as<int32_t>());
    }
    if (vm.count(kSmoothMu)) {
        param.setSmoothMu(vm[kSmoothMu].as<double>());
    }
    if (vm.count(kSeqError)) {
        param.setSeqError(vm[kSeqError].as<double>());
    }
    if (vm.count(kEmitError)) {
        param.setEmitError(vm[kEmitError].as<double>());
    }
    if (vm.count(kTranNe)) {
        param.setTranNe(vm[kTranNe].as<double>());
    }
    if (vm.count(kVerbose)) {
        param.setVerbose(vm[kVerbose].as<uint8_t>());
    }
}

int main(int argc, char** argv) {

  printVersionInformation();
  srand(time(NULL));

  po::variables_map vm;
  try {
    
    po::options_description desc("Haplotype Phasing Tool (Harsh)");
    
    desc.add_options()
    (kSeqFile, po::value<std::string>(), "sequencing read file")
    (kMapFile, po::value<std::string>(), "SNP map file")
    (kOutFile, po::value<std::string>(), "output file")
    (kRefFile, po::value<std::string>(), "reference haplotype file")
    (kMaxIter, po::value<int32_t>(), " number of sampling iterations (default 1e3)")
    (kSmoothMu, po::value<double>(), "smoothing parameter for sampling (default 1)")
    (kSeqError, po::value<double>(), "sequencing error rate (default 0.01)")
    (kEmitError, po::value<double>(), "mismatch error rate between reference and donar haplotype (default 2e-3)")
    (kTranNe, po::value<double>(), "population constant for computing transitions (default 1e-3)")
    (kVerbose, po::value<uint8_t>()->default_value(1), "verbose level (default 1)l")
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (argc == 1) {
      std::cout << desc << std::endl;
      return 1;
    }
  } catch (const std::exception &exception) {
    // Catch the initialization error and bail immediately
    std::cerr << "Error: " << exception.what() << std::endl;
  }
  
  Parameter param;
  parseInputParameter(vm, param);
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

