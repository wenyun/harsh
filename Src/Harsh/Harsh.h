#pragma once

#include <cstdint>
#include <vector>

using std::vector;
using std::pair;

class Haplotype;
class Parameter;
class Read;
class Reference;

void exitWithHelp();
void printVersionInformation();
void printInputInformation(const vector<Read>& read,
                           const Reference& reference,
                           const Parameter& param);

void printInputInformation(const vector<Read>& read,
                           const Parameter& param);

void parseInputParameter(int argc,
                         char** argv,
                         Parameter& param);

void assemblyGibbs(Haplotype& haplotype,
                   vector<Read>& read,
                   const Parameter& param);

int assemblyMismatch(Haplotype& haplotype);

void inferenceGibbs(Haplotype& haplotype1,
                    Haplotype& haplotype2,
                    vector<Read>& read,
                    const Reference& reference,
                    const Parameter& param);

double inferenceLikelihood(Haplotype& haplotype1,
                           Haplotype& haplotype2,
                           vector<Read>& read,
                           const Reference& reference,
                           const Parameter& param);

void warmStart(Haplotype& haplotype1,
               Haplotype& haplotype2,
               const vector<pair<uint16_t, uint16_t> >& count,
               const Reference& reference,
               const Parameter& param);

void plusEmission(vector<double>& v,
                  const int i,
                  const pair<uint16_t, uint16_t>& count,
                  const Reference& reference,
                  const Parameter& param);

void plusTransition(vector<double>& tmp,
                    const double extra,
                    const int j,
                    const int k,
                    const int hapNo);

int findMax(const vector<double>& score, 
            const int i, 
            const int hapNo, 
            double* maxScore);
