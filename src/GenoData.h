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

#ifndef LMMNET_DATAUTILS_H
#define LMMNET_DATAUTILS_H

#include <vector>
#include <string>
#include <map>
#include <boost/utility.hpp>

#include "TypeDef.h"
#include "IOUtils.h"
#include "InfoStructure.h"
#include "GenotypeBasis.h"

namespace LMMNET {
class GenoData : public GenoBasis {

 private:
  uchar *genotypes; // M x Nstride / 4 genotype data

  std::map<std::string, uint64> snpID_position; // to store the snpID and position

  std::vector<SnpInfo> processSnps(std::vector<uint64> &Mfiles,
                                   const std::vector<std::string> &bimFiles,
                                   const std::vector<std::string> &removeSNPsFiles);

 public:
  GenoData(const std::string &_famFile, const std::vector<std::string> &_bimFiles,
           const std::vector<std::string> &_bedFiles,
           const std::vector<std::string> &_removeSNPsFiles,
           const std::vector<std::string> &_removeIndivsFiles, double _maxMissingPerSnp,
           double _maxMissingPerIndiv);

  GenoData(const GenoData &) = delete; // disable copy constructor
  ~GenoData();

  void decodeSnpsVector(double out[], const double maskIndivs[], uint64 m,
                        const double map0129[4], double (*workTable)[4]) const;
  const std::map<std::string, uint64> &getsnpRef() const;
};
}
#endif //LMMNET_DATAUTILS_H
