// Mixture Model Net
// Copyright (C) 2019  Shunkang Zhang
// Copyright (C) 2014-2017  Harvard University
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef LMMNET_LMMCPU_H
#define LMMNET_LMMCPU_H

#include "GenoData.h"
#include "CovarBasis.h"

#include <mpi.h>

namespace LMMNET {
class LMMCPU {

 private:
  const GenoData &genoData; // class for genotype data
  const CovarBasis<GenoData> &covarBasis; // class for covariate data

  const double *maskIndivs; // final individual mask

  uint64 M, Npad;
  uint64 Nused;

  uint64 M_MPI; // total snps for MPI version

  int numChrom; // the number of chrom used in model for leave out strategy
  int numCalibSnps;
  int snpsPerBlock; // split the matrix into block to do computation
  uint64 estIteration; // B in the initial algorithm

  uchar *projMaskSnps;

  double sigma2g; // variance of signal
  double sigma2e; // variance of noise

  double *L2normXcols; //
  double (*snpLookupTable)[4]; // table for looking up raw genotype data

  bool useExactTrace; // whether compute the exact trace

  double *conjugateResultFixEff;

  std::vector<double> projGenoStd; // project genotype standard deviation
  std::vector<double> VinvyNorm; // l2 norm of Conjugate gradient result
  std::vector<double> fixEffect;
  std::vector<double> posteriorMean;

  double subIntercept;

  std::vector<std::pair<double, double> > Meanstd; // genotype mean and std

  std::vector<int> chromID;
  std::map<int, uint64> numSnpsPerChrom;
  std::vector<uint64> chromStartPoint;

  std::string imputeMethod;

  const std::string outputFile;

  // statistics used to calculate association test
  double calibrationFactor;
  std::vector<double> z;
  std::vector<double> zsquare;

  // initialize the lmmcpu object and normalize snps
  void initialize();
  uchar normalizeSnps(uint64 m, double *snpVector);
  void invnormalizeSnps();

  // decode snps
  inline void buildMaskedSnpCovCompVec(double snpCovCompVec[], uint64 m, double (*work)[4]) const {
    genoData.decodeSnpsVector(snpCovCompVec, maskIndivs, m, snpLookupTable[m], work);
  }

  // compute the heritability
  double calyVKVy(const double *projectPheno) const;
  void multXXT(double *out, const double *vec) const;
  void computeXXT(double *out) const;
  void multXXTTrace(double *out, const double *vec) const; // moment estimate for the trace
  void calTraceMoM(double &kv, double &kvkv) const;
  void calTraceExact(double &kv, double &kvkv) const;
  double calStandError(const double *projectPheno, double kv, double kvkv) const;

  // compute the association test statistics
  void makeConjugateMask(uchar *mask); // make masks for the leave out strategy
  void multXTmatrix(double *out, const double *matrix, int col) const;
  void multXmatrix(double *out, const double *matrix, int col) const;
  void multXXTConjugate(double *out, const double *matrix, const uchar *mask);
  void calConjugateBatch(double *VinvChromy, const double *inputMatrix); // compute batch conjguate gradient
  void calBeta(const double *VinvChromy); // compute beta with leave out strategy
  std::vector<uint64> selectProSnps(int numCalibSnps);

  // compute the estimate of fixed effect based on training data
  void multXXTConjugateWithoutMask(double *out, const double *matrix, unsigned int batchsize);
  void calConjugateWithoutMask(double *Viny, const double *inputVec, int batchsize);

  // overload blas interface
  double dotVec(const double *vec1, const double *vec2, uint64 N) const;
  void scalVec(double *vec, double alpha, uint64 elem) const;
  void compInverse(double *matrix,
                   const unsigned int row) const; // Note that this function will overwrite original matrix directly
  // this function imputes the missing value by using mode instead of mean
  // Todo: we should also provide the mean value to impute the data
  void buildPredictSnpsLookupTable(uint64 m,
                                   uint64 numPredict,
                                   double *snpVector,
                                   double (*predictionSnpLookupTable)[4]) const;
  void predictRandomEff(uint64 numPredict,
                        double *predictMaskIndivs,
                        double *randomEff,
                        double (*predictionSnpLookupTable)[4],
                        const GenoData &predictData) const;
  void predictFixEff(uint64 numPredict, double *fixEff, const double *predictCovarMatrix) const;

  // API for MPI version
  void multXXTConjugateWithoutMask_MPI(double* out, const double* matrix, unsigned int batchsize);
  void calConjugateWithoutMask_MPI(double *Viny, const double *inputVec, int batchsize);

 public:
  LMMCPU(const GenoData &_genoData,
         const CovarBasis<GenoData> &_covarBasis,
         const double _maskIndivs[],
         int _snpsPerBlock,
         uint64 _estIteration,
         int _numChrom,
         int _numCalibSnps,
         bool _useExactTrace,
         const std::string _imputeMethod,
         const std::string _outputFile);
  ~LMMCPU();

  void calHeritability(const double *projectPheno);
  void calCalibrationFactor(const double *projectPheno, bool sampleCalibrationFactor);
  void computeStatistics(std::string &outputFile) const;
  void estimateFixEff(const double *projectPheno, bool useApproximate);
  void computePosteriorMean(const double* pheno, bool useApproximate);
  void predict(double *output, const GenoData &predictData, const CovarBasis<GenoData> &predictCov) const;

  // API for MPI version
  void calHeritability_MPI(const double *projectPheno);
  void estimateFixEff_MPI(const double* Pheno, bool useApproximate);
  void computePosteriorMean_MPI(const double* pheno, bool useApproximate);
};
}
#endif //LMMNET_LMMCPU_H
