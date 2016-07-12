#pragma once

#include <vector>

using std::vector;

double logsum(const double logX, const double logY);

double logsumArray(const vector<double>& logArray,
                   const int start,
                   const int length);

int mnrnd(const vector<double>& prob);

double perm(int n);

double nChoosek(int n, int k);

void nSamplek(vector<int>& samples, int n, int k);
