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

  // split bed files
  int numbed_perprocess = static_cast<int>(params.bedFiles.size() + total_process - 1) / total_process;
  if (numbed_perprocess == 0) {
    cerr << "Error: the number of process " << total_process << " is larger than number of bedFiles "
    << params.bedFiles.size() << endl;
    cerr << "You should use less process" << endl;
    exit(1);
  }

  auto start_id = id * numbed_perprocess;
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

  lmmcpu.calHeritability_MPI(phenodbl[0].data());
  MPI_Finalize();
  if (id == 0) {
    cout << "Total elapsed time for analysis = " << (timer.get_time() - start_time) << " sec"
         << endl;
  }
  return 0;
}