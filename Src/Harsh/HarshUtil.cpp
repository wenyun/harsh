#include "HarshUtil.h"

#include "HarshIO.h"

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <vector>

using std::vector;

double logsum(double logX, double logY) {
  if (logX > logY) {
    return logX + log(1.0 + exp(logY - logX));
  } else {
    return logY + log(1.0 + exp(logX - logY));
  }
}

double logsumArray(const vector<double>& logArray, 
                   const int start,
                   const int length) {

  double logSum = -DBL_MAX;
  for (int i = start; i < start + length; i++) {
    logSum = logsum(logSum, logArray[i]);
  }
  return logSum;
}

int mnrnd(const vector<double>& prob) {
  int index = -1;
  while (index < 0) {
    double r = rand() / (double) RAND_MAX;
    double sum = 0.0;
    for (size_t i = 0; i < prob.size(); i++) {
      sum += prob[i];
      if (sum >= r) {
        index = i + 1;
        break;
      }
    }
  }

  return index;
}

double perm(int n) {
  int i;

  double res = 1;
  for (i = n; i > 1; i--) {
    res =  res * i;
  }
  return res;
}

double nChoosek(int n, int k) {
  if (n < 0 || k < 0) {
    exit(1);
  }
  double np = perm(n);
  double kp = perm(k);
  double nkp = perm(n-k);
  
  return np / (kp * nkp);
}

void nSamplek(vector<int>& samples, int n, int k) {

  samples.resize(k);
  int i, j;
  for (i = 0, j = 0; i < n && j < k; ++i) {
    int ri = n - i;
    int rj = k - j;
    if (rand() % ri < rj) {
      samples[j++] = i;
    }
  }
  
  if (j != k) {
    harshErrorExit("Knuth algorithm wrong\n");
  }
}
