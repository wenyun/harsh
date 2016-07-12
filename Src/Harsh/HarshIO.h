#pragma once

#include <string>
#include <vector>

using std::string;
using std::vector;

class Read;
class Parameter;
class Haplotype;

void readSeqFile(const char* seqFile, 
                 vector<Read>& read,
                 const Parameter& param);
void outputAssembly(const Haplotype& haplotype,
                    const Parameter& param);
void outputInference(const Haplotype& haplotype1,
                     const Haplotype& haplotype2,
                     const Parameter& param);
bool isEmpty(const string& line);
void harshErrorExit(const char* msg);
void harshWarning(const char* msg);
void harshMessage(const char* msg, 
                  const int level,
                  const Parameter& param);
