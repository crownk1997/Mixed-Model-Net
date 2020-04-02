//
// Created by zsk on 2020/3/12.
//

#include "GenoData.h"
#include "PhenoBasis.h"
#include "CovarBasis.h"
#include "Parameters.h"
#include "LMMCPU.h"
#include "Timer.h"

#include <mpi.h>
#include "omp.h"
#include <mkl.h>

#include <iostream>

using namespace std;
using namespace LMMNET;

int main(int argc, char** argv) {

  Timer timer;
  double start_time = timer.get_time();

  MPI_Init(&argc, &argv);

  int total_process;
  int id;
  MPI_Comm_size(MPI_COMM_WORLD, &total_process);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if (id == 0) {
    cout << "========== Testing MPI properties for MMNET ==========" << endl;
  }

  // the rank 0 output parameters
  if (id == 0) {
    cout << "The input arguments are :" << "\n";
    for (int num = 1; num < argc; num++) {
      cout << argv[num++] << " " << argv[num] << "\n";
    }
    cout << endl;
  }

  Params params;
  if (!params.processCommandLineArgs(argc, argv)) {
    cerr << "Aborting due to error processing command line arguments" << endl;
    cerr << "For list of arguments, run with -h (--help) option" << endl;
    exit(1);
  }

  if (id == 0) {
    cout << endl << "Setting number of threads to " << params.numThreads << endl;
  }

  omp_set_num_threads(params.numThreads);
  mkl_set_num_threads(params.numThreads);

  if (id == 0) {
    cout << "fam: " << params.famFile << endl;
    cout << "bim(s): ";
    for (uint i = 0; i < params.bimFiles.size(); i++) cout << params.bimFiles[i] << endl;
    cout << "bed(s): ";
    for (uint i = 0; i < params.bedFiles.size(); i++) cout << params.bedFiles[i] << endl;
  }

  /********** SET UP DATA **********/
  if (id == 0) {
    cout << endl << "***** Reading genotype data *****" << endl << endl;
  }

  // According to total number of processes to split input file
  int total_files = params.bedFiles.size();
  int numbed_perprocess;
  int start_id;
  if (total_files < total_process) {
    cerr << "Error: the number of process " << total_process << " is larger than number of bedFiles "
    << total_files << endl;
    cerr << "You should use less process" << endl;
    exit(1);
  }

  int remainder = total_files % total_process;
  if (remainder != 0) {
    cout << "Warning: input files cannot be seperated evenly " << endl;
    cout << "There might be imbalance workload among different processes" << endl;
    // compute the remainder and split task evenly
//    int zp = total_files - remainder;
    int pp = total_files / total_process;

    if (id < remainder) {
      numbed_perprocess = pp + 1;
      start_id = id * numbed_perprocess;
    } else {
      numbed_perprocess = pp;
      start_id = remainder * (pp + 1) + (id - remainder) * numbed_perprocess;
    }
  } else {
    // task can be divided evenly
    numbed_perprocess = total_files / total_process;
    start_id = id * numbed_perprocess;
  }

  auto end_id = start_id + numbed_perprocess;
  const std::vector<std::string> bedFiles_PerP(&params.bedFiles[start_id], &params.bedFiles[end_id]);
  const std::vector<std::string> bimFiles_PerP(&params.bimFiles[start_id], &params.bimFiles[end_id]);

  
  GenoData genoData(params.famFile, bimFiles_PerP, bedFiles_PerP, params.removeSnpsFiles,
                    params.removeIndivsFiles, params.maxMissingPerSnp, params.maxMissingPerIndiv);
  const vector<SnpInfo> &snps = genoData.getSnpInfo();

  cout << "Time for setting up dataset and read data " << timer.update_time() << " sec" << endl;

  if ((int) snps.size() > params.maxModelSnps) {
    cerr << "ERROR: Number of SNPs exceeds maxModelSnps = " << params.maxModelSnps << endl;
    exit(1);
  }

  cout << endl << "***** Reading phenotype and covariate data *****" << endl << endl;

  vector<double> maskIndivs(genoData.getNpad());
  genoData.writeMaskIndivs(maskIndivs.data());

  // read phenotype data from files if provided
  PhenoBasis<GenoData> phenobasis(params.phenoFile, genoData, params.phenoCols, maskIndivs, params.phenoUseFam);
  phenobasis.padPhenodbl(); // pad the phenotype data to have the same individuals as genotype data

  vector<vector<double> > phenodbl = phenobasis.getPhenodbl();

  // read covariate data from files if provided
  CovarBasis<GenoData> covarBasis(params.covarFile, genoData, params.covarCols, maskIndivs);
  covarBasis.writeMask(maskIndivs);
  // update the phenotype value with new mask
  for (uint64 n = 0; n < phenodbl[0].size(); n++) {
    phenodbl[0][n] *= maskIndivs[n];
  }
  covarBasis.projectCovarsVec(phenodbl[0].data());

  cout << "Time for reading and processing covariate and phenotype data " << timer.update_time() << " esc" << endl;

  cout << endl << "***** Initializing lmmnet object and normlizing snps *****" << endl << endl;

  LMMCPU lmmcpu(genoData, covarBasis, &maskIndivs[0], params.snpsPerBlock, params.estIterationTrace, params.numChrom,
                params.numCalibSnps, params.useExactTrace, params.imputeMethod, params.outputFile);

  cout << "Time for initializing lmmnet object and normalizing snps is " << timer.update_time() << " sec" << endl;

  cout << endl << "***** Computing heritability *****" << endl << endl;

  lmmcpu.calHeritability_MPI(phenodbl[0].data());

  cout << "Time for computing heritability is " << timer.update_time() << " sec" << endl;

  if (params.prediction) {
    cout << endl << "***** Compute the estimate for fix effect w *****" << endl << endl;

    vector<vector<double> > phenodblorigin = phenobasis.getPhenodbl(); // here we use the original pheno data
    lmmcpu.estimateFixEff_MPI(phenodblorigin[0].data(), params.useApproFixEffect);

    cout << "Timer for estimating fix effect " << timer.update_time() << " esc" << endl;

    cout << endl << "***** Compute the posterior mean *****" << endl << endl;
    lmmcpu.computePosteriorMean(phenodblorigin[0].data(), params.useApproFixEffect);

    cout << "Timer for computing posterior mean " << timer.update_time() << " esc" << endl;

    cout << endl << "***** Read the prediction data *****" << endl << endl;
    GenoData singlepredictData(params.predfamFile, params.predbimFiles, params.predbedFiles, params.predremoveSnpsFiles,
                               params.predremoveIndivsFiles, params.maxMissingPerSnp, params.maxMissingPerIndiv);

    const vector<SnpInfo> &predictSnps = singlepredictData.getSnpInfo();

    // the number of prediction SNPs cannot be too large
    if (predictSnps.size() > params.maxModelSnps) {
      cerr << "ERROR: Number of prediction SNPs exceeds maxModelSnps = " << params.maxModelSnps << endl;
      exit(1);
    }

    // the number of prediction snps should be equal to training snps
    if (predictSnps.size() != snps.size()) {
      cerr << "Error: Number of prediction SNPs should be equal the number of training SNPs " << endl;
      exit(1);
    }

    // read predict covariate data
    vector<double> singPredmaskIndivs(singlepredictData.getNpad());
    singlepredictData.writeMaskIndivs(singPredmaskIndivs.data());
    CovarBasis<GenoData> predictCov(params.predcovarFile, singlepredictData, params.precovarCols, singPredmaskIndivs);

    cout << "Time for setting up dataset and read prediction data " << timer.update_time() << " sec" << endl;

    // predict the new phenotype based on the posterior mean

    double *predictOutput = ALIGN_ALLOCATE_DOUBLES(singlepredictData.getNpad());
    lmmcpu.predict(predictOutput, singlepredictData, predictCov);

    cout << "Timer for predict new data " << timer.update_time() << " esc" << endl;

  }
  MPI_Finalize();
  if (id == 0) {
    cout << "Total elapsed time for analysis = " << (timer.get_time() - start_time) << " sec"
         << endl;
  }
  return 0;
}