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

#ifndef LMMNET__GENOTYPEBASIS_H_
#define LMMNET__GENOTYPEBASIS_H_

#include <vector>
#include <string>
#include <map>
//#include <boost/utility.hpp>

#include "TypeDef.h"
#include "IOUtils.h"
#include "InfoStructure.h"

namespace LMMNET {
class GenoBasis {
 public:
  static const uint64 IND_MISSING;

 protected:
  uint64 Mbed, Nbed; // PLINK data dimensions
  std::vector<int> bedIndivToRemoveIndex; // index of removing

  uint64 M, N, Nstride; // PLINK data dimensions after data processing
  uint64 Nused;

  uchar *maskSnps;
  double *maskIndivs; // mask of indivs failed QC
  uint64 numIndivsQC; // number of indivs remaining after QC

  std::map<std::string, uint64> FID_IID_to_ind; // map FID+IID to ind

  double maxMissingPerIndiv; // max missing rate per individual
  double maxMissingPerSnp; // max missing rate per snp

  // path for bim and bed files
  std::vector<std::string> bimFiles;
  std::vector<std::string> bedFiles;

  // container for indivs and snps information
  std::vector<SnpInfo> snps;
  std::vector<IndivInfo> indivs;

  void storeBedLine(uchar bedLineOut[], const uchar genoLine[]);
  void processIndivs(const std::string &famFile, const std::vector<std::string> &removeIndivsFiles);

 public:
  GenoBasis(const std::string &famFile, const std::vector<std::string> &bimFiles,
            const std::vector<std::string> &bedFiles,
            const std::vector<std::string> &removeSNPsFiles,
            const std::vector<std::string> &removeIndivsFiles, double _maxMissingPerSnp,
            double _maxMissingPerIndiv);

  GenoBasis(const GenoBasis &) = delete; // disable copy constructor
  ~GenoBasis();

  static std::vector<SnpInfo> readBimFile(const std::string &bimFile);

  void readBedLine(uchar genoLine[], uchar bedLineIn[], FileUtils::SafeIfstream &fin) const;
  double computeAlleleFreq(const uchar genoLine[], const double subMaskIndivs[]) const;
  double computeMAF(const uchar genoLine[], const double subMaskIndivs[]) const;
  double computeSnpMissing(const uchar genoLine[], const double subMaskIndivs[]) const;

  void buildLookupTable(double (*workTable)[4], const double lookupBedTable[4]) const;

  // interface of get and set
  uint64 getIndivInd(std::string &FID, std::string &IID) const;
  const std::vector<SnpInfo> &getSnpInfo() const;
  void updateNused();

  uint64 getNpad() const;
  uint64 getM() const;
  uint64 getNused() const;
  std::vector<double> getFamPhenos() const;
  // output individual and snps masks to file
  void writeMaskSnps(uchar out[]) const;
  void writeMaskIndivs(double out[]) const;
};
}

#endif //LMMNET__GENOTYPEBASIS_H_
