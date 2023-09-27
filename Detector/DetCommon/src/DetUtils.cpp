#include "DetCommon/DetUtils.h"

// DD4hep
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4VolumeManager.h"

// Geant
#include "G4NavigationHistory.hh"

// ROOT
#include "TGeoBBox.h"

#ifdef HAVE_GEANT4_UNITS
#define MM_2_CM 1.0
#else
#define MM_2_CM 0.1
#endif

// GM for debug
#include <iostream>

#include <unordered_set>

namespace det {
namespace utils {
dd4hep::xml::Component getNodeByStrAttr(const dd4hep::xml::Handle_t& mother, const std::string& nodeName,
                                        const std::string& attrName, const std::string& attrValue) {
  for (dd4hep::xml::Collection_t xCompColl(mother, nodeName.c_str()); nullptr != xCompColl; ++xCompColl) {
    if (xCompColl.attr<std::string>(attrName.c_str()) == attrValue) {
      return static_cast<dd4hep::xml::Component>(xCompColl);
    }
  }
  // in case there was no xml daughter with matching name
  return dd4hep::xml::Component(nullptr);
}

double getAttrValueWithFallback(const dd4hep::xml::Component& node, const std::string& attrName,
                                const double& defaultValue) {
  if (node.hasAttr(_Unicode(attrName.c_str()))) {
    return node.attr<double>(attrName.c_str());
  } else {
    return defaultValue;
  }
}

uint64_t cellID(const dd4hep::Segmentation& aSeg, const G4Step& aStep, bool aPreStepPoint) {
  dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
  dd4hep::VolumeID volID = volMgr.volumeID(aStep.GetPreStepPoint()->GetTouchable());
  if (aSeg.isValid()) {
    G4ThreeVector global;
    if (aPreStepPoint) {
      global = aStep.GetPreStepPoint()->GetPosition();
    } else {
      global = 0.5 * (aStep.GetPreStepPoint()->GetPosition() + aStep.GetPostStepPoint()->GetPosition());
    }
    G4ThreeVector local =
        aStep.GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(global);
    dd4hep::Position loc(local.x() * MM_2_CM, local.y() * MM_2_CM, local.z() * MM_2_CM);
    dd4hep::Position glob(global.x() * MM_2_CM, global.y() * MM_2_CM, global.z() * MM_2_CM);
    dd4hep::VolumeID cID = aSeg.cellID(loc, glob, volID);
    return cID;
  }
  return volID;
}

std::vector<std::vector<uint>> combinations(int N, int K) {
  std::vector<std::vector<uint>> indexes;
  std::string bitmask(K, 1);  // K leading 1's
  bitmask.resize(N, 0);       // N-K trailing 0's
  // permute bitmask
  do {
    std::vector<uint> tmp;
    for (int i = 0; i < N; ++i) {  // [0..N-1] integers
      if (bitmask[i]) {
        tmp.push_back(i);
      }
    }
    indexes.push_back(tmp);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return indexes;
}

std::vector<std::vector<int>> permutations(int K) {
  std::vector<std::vector<int>> indexes;
  int N = pow(2, K);  // number of permutations with repetition of 2 numbers (-1,1)
  for (int i = 0; i < N; i++) {
    // permutation = binary representation of i
    std::vector<int> tmp;
    tmp.assign(K, 0);
    uint res = i;
    // dec -> bin
    for (int j = 0; j < K; j++) {
      tmp[K - 1 - j] = -1 + 2 * (res % 2);
      res = floor(res / 2);
    }
    indexes.push_back(tmp);
  }
  return indexes;
}

// use it for field module/phi
int cyclicNeighbour(int aCyclicId, std::pair<int, int> aFieldExtremes) {
  int nBins = aFieldExtremes.second - aFieldExtremes.first + 1;
  if (aCyclicId < aFieldExtremes.first) {
    return  (aFieldExtremes.second + 1 - ((aFieldExtremes.first - aCyclicId) % nBins) );
  } else if (aCyclicId > aFieldExtremes.second) {
    return ( ((aCyclicId - aFieldExtremes.first) % nBins) + aFieldExtremes.first) ;
  }
  return aCyclicId;
}

std::vector<uint64_t> neighbours(const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::pair<int, int>>& aFieldExtremes, uint64_t aCellId,
                                 const std::vector<bool>& aFieldCyclic, bool aDiagonal) {
  std::vector<uint64_t> neighbours;
  dd4hep::DDSegmentation::CellID cID = aCellId;
  for (uint itField = 0; itField < aFieldNames.size(); itField++) {
    const auto& field = aFieldNames[itField];
    // note: get(..) returns a FieldID i.e. a int64_t
    // while CellID (used previously) is a uint64_t
    // similarly, the second argument of set(..)
    // is signed (FieldID) rather than unsigned (CellID)
    int id = aDecoder.get(cID,field);
    if (aFieldCyclic[itField]) {
      aDecoder[field].set(cID, cyclicNeighbour(id - 1, aFieldExtremes[itField]));
      neighbours.emplace_back(cID);
      aDecoder[field].set(cID, cyclicNeighbour(id + 1, aFieldExtremes[itField]));
      neighbours.emplace_back(cID);
    } else {
      if (id > aFieldExtremes[itField].first) {
        aDecoder[field].set(cID, id - 1);
        neighbours.emplace_back(cID);
      }
      if (id < aFieldExtremes[itField].second) {
        aDecoder[field].set(cID, id + 1);
        neighbours.emplace_back(cID);
      }
    }
    aDecoder[field].set(cID, id);
  }
  if (aDiagonal) {
    std::vector<int> fieldIds;  // initial IDs
    fieldIds.assign(aFieldNames.size(), 0);
    // for each field get current Id
    for (uint iField = 0; iField < aFieldNames.size(); iField++) {
      const auto& field = aFieldNames[iField];
      fieldIds[iField] = aDecoder.get(cID, field);
    }
    for (uint iLength = aFieldNames.size(); iLength > 1; iLength--) {
      // get all combinations for a given length
      const auto& indexes = combinations(aFieldNames.size(), iLength);
      for (uint iComb = 0; iComb < indexes.size(); iComb++) {
        // for current combination get all permutations of +- 1 operation on IDs
        const auto& calculation = permutations(iLength);
        // do the increase/decrease of bitfield
        for (uint iCalc = 0; iCalc < calculation.size(); iCalc++) {
          // set new Ids for each field combination
          bool add = true;
          for (uint iField = 0; iField < indexes[iComb].size(); iField++) {
            if (aFieldCyclic[indexes[iComb][iField]]) {
              aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, cyclicNeighbour(fieldIds[indexes[iComb][iField]] + calculation[iCalc][iField],
                                                                                     aFieldExtremes[indexes[iComb][iField]]) );
            } else if ((calculation[iCalc][iField] > 0 &&
                        fieldIds[indexes[iComb][iField]] < aFieldExtremes[indexes[iComb][iField]].second) ||
                       (calculation[iCalc][iField] < 0 &&
                        fieldIds[indexes[iComb][iField]] > aFieldExtremes[indexes[iComb][iField]].first)) {
              aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, fieldIds[indexes[iComb][iField]] + calculation[iCalc][iField]);
            } else {
              add = false;
            }
          }
          // add new cellId to neighbours (unless it's beyond extrema)
          if (add) {
            neighbours.emplace_back(cID);
          }
          // reset ids
          for (uint iField = 0; iField < indexes[iComb].size(); iField++) {
            aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, fieldIds[indexes[iComb][iField]]);
          }
        }
      }
    }
  }
  return neighbours;
}

// use it for module-theta merged readout (FCCSWGridModuleThetaMerged)
std::vector<uint64_t> neighbours_ModuleThetaMerged(const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged& aSeg,
                                 const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                 const std::vector<std::string>& aFieldNames,
                                 const std::vector<std::pair<int, int>>& aFieldExtremes, uint64_t aCellId,
                                 const std::vector<bool>& aFieldCyclic, bool aDiagonal) {
  std::vector<uint64_t> neighbours;
  int nLayer = aDecoder.get(aCellId, aSeg.fieldNameLayer());
  dd4hep::DDSegmentation::CellID cID = aCellId;

  for (uint itField = 0; itField < aFieldNames.size(); itField++) {
    // retrieve name of field and corresponding value in cellID
    const auto& field = aFieldNames[itField];
    //dd4hep::DDSegmentation::CellID id = aDecoder.get(aCellId, field);
    int id = aDecoder.get(aCellId, field);

    if (field == aSeg.fieldNameLayer()) {
      // for neighbours across different layers, we have to take into
      // account that merging along module and/or theta could be different
      // so one cell in layer N could be neighbour to several in layer N+-1
      //dd4hep::DDSegmentation::CellID module_id = aDecoder.get(aCellId, aSeg.fieldNameModule());
      //dd4hep::DDSegmentation::CellID theta_id = aDecoder.get(aCellId, aSeg.fieldNameTheta());
      int module_id = aDecoder.get(aCellId, aSeg.fieldNameModule());
      int theta_id = aDecoder.get(aCellId, aSeg.fieldNameTheta());

      /*
      // for inner layer (N-1)
      if (id > aFieldExtremes[itField].first) {
        aDecoder.set(cID, field, id - 1);
        if (aSeg.mergedModules(nLayer) == aSeg.mergedModules(nLayer-1)) {
	  // if the numbers of merged cells across 2 layers are the same we just do +-1
          neighbours.emplace_back(cID);
        } else {
	  // otherwise, we need to do some math
	  // note: even if number of merged cells in layer-1 is larger, a cell
	  // in layer could neighbour more than one cell in layer-1 if the merged
	  // cells are not aligned, for example if cells are grouped by 3 in a layer
	  // and by 4 in the next one, cell 435 in the former (which groups together
	  // 435-436-437) will be neighbour to cells 432 and 436 of the latter
	  // this might introduce duplicates, we will remove them later
	  // another issue is that it could add spurious cells beyond the maximum module number
	  // to prevent this we would need to know the max module number in layer -1
	  // which would require modifying this function passing the extrema for all layers
	  // instead of the extrema only for a certain layer
	  // this border effect is also present in the original method..
          for (int i=0; i < aSeg.mergedModules(nLayer); i++) {
            aDecoder.set(cID, aSeg.fieldNameModule(), (module_id + i) - ((module_id + i) % aSeg.mergedModules(nLayer-1)));
            neighbours.emplace_back(cID);
	  }
        }
	// reset module ID
	aDecoder.set(cID, aSeg.fieldNameModule(), module_id);

        // do the same for theta
	// this might introduce some duplicates
        if (aSeg.mergedThetaCells(nLayer) == aSeg.mergedThetaCells(nLayer-1)) {
	  // if the numbers of merged cells across 2 layers are the same we just do +-1
          neighbours.emplace_back(cID);
        } else {
          for (int i=0; i < aSeg.mergedThetaCells(nLayer); i++) {
            aDecoder.set(cID, aSeg.fieldNameTheta(), (theta_id + i) - ((theta_id + i) % aSeg.mergedThetaCells(nLayer-1)));
            neighbours.emplace_back(cID);
	  }
        }
	// reset theta
	aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id);
      }
      // for outer layer (N+1)
      if (id < aFieldExtremes[itField].second) {
        aDecoder.set(cID, field, id + 1);
        if (aSeg.mergedModules(nLayer) == aSeg.mergedModules(nLayer+1)) {
          neighbours.emplace_back(cID);
	  // note: in theta, there can be border effects leading to cells outside of
	  // the detector volume to be added - consider the case of a theta cell at z=zmax
	  // then in outer layer the cell with same thetaID could be beyond zmax
	  // anyway these corresponds to regions outside the acceptance of the full calo
	  // so probably in clusters that get dropped because outside of fiducial region
        } else {
	  for (int i=0; i < aSeg.mergedModules(nLayer); i++) {
            aDecoder.set(cID, aSeg.fieldNameModule(), (module_id + i) - ((module_id + i) % aSeg.mergedModules(nLayer+1)));
            neighbours.emplace_back(cID);
	  }
        }
	// reset module ID
	aDecoder.set(cID, aSet.fieldNameModule(), module_id);
	// do the same for theta
	// this might introduce some duplicates
        if (aSeg.mergedThetaCells(nLayer) == aSeg.mergedThetaCells(nLayer+1)) {
          neighbours.emplace_back(cID);
        } else {
          for (int i=0; i < aSeg.mergedThetaCells(nLayer); i++) {
            aDecoder.set(cID, aSeg.fieldNameTheta(), (theta_id + i) - ((theta_id + i) % aSeg.mergedThetaCells(nLayer+1)));
            neighbours.emplace_back(cID);
	  }
        }
	// reset theta
	aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id);
      }
    }
    */
      for (int deltaLayer = -1; deltaLayer<2; deltaLayer+=2) {
	// no neighbours in layer N-1 for innermost layer
	if (id == aFieldExtremes[itField].first && deltaLayer<0) continue;
	// and in layer N+1 for outermost layer
	if (id == aFieldExtremes[itField].second && deltaLayer>0) continue;

	// set layer field of neighbour cell
        aDecoder.set(cID, field, id + deltaLayer);

	// find the neighbour(s) in module and theta
	// if the numbers of module (theta) merged cells across 2 layers are the
	// same then we just take the same module (theta) ID
	// otherwise, we need to do some math to account for the different mergings
	// note: even if number of merged cells in layer-1 is larger, a cell
	// in layer could neighbour more than one cell in layer-1 if the merged
	// cells are not aligned, for example if cells are grouped by 3 in a layer
	// and by 4 in the next one, cell 435 in the former (which groups together
	// 435-436-437) will be neighbour to cells 432 and 436 of the latter
	// this might introduce duplicates, we will remove them later
	// another issue is that it could add spurious cells beyond the maximum module number
	// to prevent this we would need to know the max module number in layer -1
	// which would require modifying this function passing the extrema for all layers
	// instead of the extrema only for a certain layer
	// this border effect is also present in the original method..
	std::vector<int> neighbourModules;
        if (aSeg.mergedModules(nLayer) == aSeg.mergedModules(nLayer+deltaLayer)) {
          neighbourModules.emplace_back(module_id);
        }
	else {
          for (int i=0; i < aSeg.mergedModules(nLayer); i++) {
            neighbourModules.emplace_back((module_id + i) - ((module_id + i) % aSeg.mergedModules(nLayer+deltaLayer)));
	  }
        }
	std::vector<int> neighbourThetaCells;
        if (aSeg.mergedThetaCells(nLayer) == aSeg.mergedThetaCells(nLayer+deltaLayer)) {
	  neighbourThetaCells.emplace_back(theta_id);
	}
	else {
	  for (int i=0; i < aSeg.mergedThetaCells(nLayer); i++) {
	    neighbourThetaCells.emplace_back((theta_id + i) - ((theta_id + i) % aSeg.mergedThetaCells(nLayer+deltaLayer)));
	  }
	}
	for (unsigned int i=0 ; i < neighbourModules.size(); i++) {
	  aDecoder.set(cID, aSeg.fieldNameModule(), neighbourModules[i]);
	  for (unsigned int j=0 ; j < neighbourThetaCells.size(); j++) {
	    aDecoder.set(cID, aSeg.fieldNameTheta(), neighbourThetaCells[j]);
	    neighbours.emplace_back(cID);
	  }
	}
      }
      // reset module and theta of cellID
      aDecoder.set(cID, aSeg.fieldNameModule(), module_id);
      aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id);
    }
    // for neighbours in module/theta direction, do +-nMergedCells instead of +-1
    else if (field == aSeg.fieldNameModule()) {
      int newid = cyclicNeighbour(id + aSeg.mergedModules(nLayer), aFieldExtremes[itField]);
      // the following line is needed in case the number of modules is not
      // a multiple of the number of merged modules
      // for instance if nModules=1545 and mergedModules =2
      // the last module, 1544, is neighbour of 0, but
      // cyclicNeighbour(1546) will return 1
      newid -= (newid % aSeg.mergedModules(nLayer));
      aDecoder[field].set(cID, newid);
      neighbours.emplace_back(cID);

      newid = cyclicNeighbour(id - aSeg.mergedModules(nLayer), aFieldExtremes[itField]);
      newid -= (newid % aSeg.mergedModules(nLayer));
      aDecoder[field].set(cID, newid);
      neighbours.emplace_back(cID);
    }
    else if (field == aSeg.fieldNameTheta()) {
      if (id-aSeg.mergedThetaCells(nLayer) >= aFieldExtremes[itField].first) {
        aDecoder[field].set(cID, id - aSeg.mergedThetaCells(nLayer));
        neighbours.emplace_back(cID);
      }
      if (id+aSeg.mergedThetaCells(nLayer) <= aFieldExtremes[itField].second) {
        aDecoder[field].set(cID, id + aSeg.mergedThetaCells(nLayer));
        neighbours.emplace_back(cID);
      }
    }
    // restore cellID
    aDecoder[field].set(cID, id);
  }

  // Diagonal case : TODO
  if (aDiagonal) {
    std::vector<int> fieldIds;  // initial IDs
    fieldIds.assign(aFieldNames.size(), 0);
    // for each field get current Id
    for (uint iField = 0; iField < aFieldNames.size(); iField++) {
      const auto& field = aFieldNames[iField];
      fieldIds[iField] = aDecoder.get(cID, field);
    }
    for (uint iLength = aFieldNames.size(); iLength > 1; iLength--) {
      // get all combinations for a given length
      const auto& indexes = combinations(aFieldNames.size(), iLength);
      for (uint iComb = 0; iComb < indexes.size(); iComb++) {
        // for current combination get all permutations of +- 1 operation on IDs
        const auto& calculation = permutations(iLength);
        // do the increase/decrease of bitfield
        for (uint iCalc = 0; iCalc < calculation.size(); iCalc++) {
          // set new Ids for each field combination
          bool add = true;
          for (uint iField = 0; iField < indexes[iComb].size(); iField++) {
            if (aFieldCyclic[indexes[iComb][iField]]) {
              aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, cyclicNeighbour(fieldIds[indexes[iComb][iField]] + calculation[iCalc][iField],
										     aFieldExtremes[indexes[iComb][iField]]) );
            } else if ((calculation[iCalc][iField] > 0 &&
                        fieldIds[indexes[iComb][iField]] < aFieldExtremes[indexes[iComb][iField]].second) ||
                       (calculation[iCalc][iField] < 0 &&
                        fieldIds[indexes[iComb][iField]] > aFieldExtremes[indexes[iComb][iField]].first)) {
              aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, fieldIds[indexes[iComb][iField]] + calculation[iCalc][iField]);
            } else {
              add = false;
            }
          }
          // add new cellId to neighbours (unless it's beyond extrema)
          if (add) {
            neighbours.emplace_back(cID);
          }
          // reset ids
          for (uint iField = 0; iField < indexes[iComb].size(); iField++) {
            aDecoder[aFieldNames[indexes[iComb][iField]]].set(cID, fieldIds[indexes[iComb][iField]]);
          }
        }
      }
    }
  }
  // remove duplicates
  std::unordered_set<uint64_t> s;
  for (uint64_t i : neighbours)
    s.insert(i);
  neighbours.assign( s.begin(), s.end() );
  return neighbours;
}

std::vector<std::pair<int, int>> bitfieldExtremes(const dd4hep::DDSegmentation::BitFieldCoder& aDecoder,
                                                  const std::vector<std::string>& aFieldNames) {
  std::vector<std::pair<int, int>> extremes;
  int width = 0;
  for (const auto& field : aFieldNames) {
    width = aDecoder[field].width();
    if (aDecoder[field].isSigned()) {
      extremes.emplace_back(std::make_pair(-(1 << (width - 1)), (1 << (width - 1)) - 1));
    } else {
      extremes.emplace_back(std::make_pair(0, (1 << width) - 1));
    }
  }
  return extremes;
}

CLHEP::Hep3Vector envelopeDimensions(uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto pvol = volMgr.lookupVolumePlacement(aVolumeId);
  auto solid = pvol.volume().solid();
  // get the envelope of the shape
  TGeoBBox* box = (dynamic_cast<TGeoBBox*>(solid.ptr()));
  // get half-widths
  return CLHEP::Hep3Vector(box->GetDX(), box->GetDY(), box->GetDZ());
}

std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXY& aSeg) {
  // // get half-widths
  auto halfSizes = envelopeDimensions(aVolumeId);
  // get segmentation cell widths
  double xCellSize = aSeg.gridSizeX();
  double yCellSize = aSeg.gridSizeY();
  // calculate number of cells, the middle cell is centred at 0 (no offset)
  uint cellsX = ceil((halfSizes.x() - xCellSize / 2.) / xCellSize) * 2 + 1;
  uint cellsY = ceil((halfSizes.y() - yCellSize / 2.) / yCellSize) * 2 + 1;
  return {cellsX, cellsY};
}

std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXYZ& aSeg) {
  // // get half-widths
  auto halfSizes = envelopeDimensions(aVolumeId);
  // get segmentation cell widths
  double xCellSize = aSeg.gridSizeX();
  double yCellSize = aSeg.gridSizeY();
  double zCellSize = aSeg.gridSizeZ();
  // calculate number of cells, the middle cell is centred at 0 (no offset)
  uint cellsX = ceil((halfSizes.x() - xCellSize / 2.) / xCellSize) * 2 + 1;
  uint cellsY = ceil((halfSizes.y() - yCellSize / 2.) / yCellSize) * 2 + 1;
  uint cellsZ = ceil((halfSizes.z() - zCellSize / 2.) / zCellSize) * 2 + 1;
  return {cellsX, cellsY, cellsZ};
}

CLHEP::Hep3Vector tubeDimensions(uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto pvol = volMgr.lookupVolumePlacement(aVolumeId);
  auto solid = pvol.volume().solid();

  // get the envelope of the shape
  TGeoTubeSeg* tube = (dynamic_cast<TGeoTubeSeg*>(solid.ptr()));
  if (tube == nullptr) {
    return CLHEP::Hep3Vector(0, 0, 0);
  }
  // get half-widths
  return CLHEP::Hep3Vector(tube->GetRmin(), tube->GetRmax(), tube->GetDZ());
}

CLHEP::Hep3Vector coneDimensions(uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto pvol = volMgr.lookupVolumePlacement(aVolumeId);
  auto solid = pvol.volume().solid();
  // get the envelope of the shape
  TGeoCone* cone = (dynamic_cast<TGeoCone*>(solid.ptr()));
  if (cone == nullptr) {
    return CLHEP::Hep3Vector(0, 0, 0);
  }
  // get half-widths
  return CLHEP::Hep3Vector(cone->GetRmin1(), cone->GetRmax1(), cone->GetDZ());
}

std::array<double, 2> tubeEtaExtremes(uint64_t aVolumeId) {
  auto sizes = tubeDimensions(aVolumeId);
  if (sizes.mag() == 0) {
    // if it is not a cylinder maybe it is a cone (same calculation for extremes)
    sizes = coneDimensions(aVolumeId);
    if (sizes.mag() == 0) {
      return {0, 0};
    }
  }
  // eta segmentation calculate maximum eta from the inner radius (no offset is taken into account)
  double maxEta = 0;
  double minEta = 0;

  double rIn = sizes.x();
  double rOut = sizes.y();
  double dZ = sizes.z();

  // check if it is a cylinder centred at z=0
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto detelement = volMgr.lookupDetElement(aVolumeId);
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  double outGlobal[3];
  double inLocal[] = {0, 0, 0};  // to get middle of the volume
  transformMatrix.LocalToMaster(inLocal, outGlobal);
  double zCenter = outGlobal[2];
  if (fabs(zCenter) < 1e-10) {
    // this assumes cylinder centred at z=0
    maxEta = CLHEP::Hep3Vector(rIn, 0, dZ).eta();
    minEta = -maxEta;
  } else {
    // GM: I think this calculation is not correct
    /*
    maxEta = CLHEP::Hep3Vector(sizes.x(), 0, std::max(sizes.z() + outGlobal[2], -sizes.z() + outGlobal[2])).eta();
    minEta = CLHEP::Hep3Vector(sizes.y(), 0, std::min(sizes.z() + outGlobal[2], -sizes.z() + outGlobal[2])).eta();
    */
    maxEta = std::max(
		      CLHEP::Hep3Vector(rIn, 0, zCenter+dZ).eta(),
		      CLHEP::Hep3Vector(rOut, 0, zCenter+dZ).eta()
		      );
    minEta = std::min(
		      CLHEP::Hep3Vector(rIn, 0, zCenter-dZ).eta(),
		      CLHEP::Hep3Vector(rOut, 0, zCenter-dZ).eta()
		      );
  }
  return {minEta, maxEta};
}

// GM: unnneded, kept for debug, will remove before merge
/*
std::array<double, 2> tubeThetaExtremes(uint64_t aVolumeId) {
  // get x-y-z sizes of cylinder
  auto sizes = tubeDimensions(aVolumeId);
  if (sizes.mag() == 0) {
    // if it is not a cylinder maybe it is a cone (same calculation for extremes)
    sizes = coneDimensions(aVolumeId);
    if (sizes.mag() == 0) {
      return {0, 0};
    }
  }
  // for some reason the code does not arrive here,
  // it looks like a cylinder is not found and thus
  // the envelope is used instead
  double rIn = sizes.x();
  double rOut = sizes.y();
  double dZ = sizes.z();
  // theta segmentation calculate maximum theta from the inner radius (no offset is taken into account)
  double maxTheta = 0;
  double minTheta = 0;
  // check if it is a cylinder centred at z=0
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto detelement = volMgr.lookupDetElement(aVolumeId);
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  double outGlobal[3];
  double inLocal[] = {0, 0, 0};  // to get middle of the volume
  transformMatrix.LocalToMaster(inLocal, outGlobal);
  double zCenter = outGlobal[2];
  if (fabs(zCenter) < 1e-6) {
    // this assumes cylinder centred at z=0
    minTheta = CLHEP::Hep3Vector(rIn, 0, dZ).theta();
    maxTheta = M_PI-minTheta;
  } else {
    // need max/min vs rOut case if origin is not inside volume of detector
    maxTheta = std::max(
			CLHEP::Hep3Vector(rIn, 0, zCenter - dZ).theta(),
			CLHEP::Hep3Vector(rOut, 0, zCenter - dZ).theta()
			);
    minTheta = std::min(
			CLHEP::Hep3Vector(rIn, 0, zCenter + dZ).theta(),
			CLHEP::Hep3Vector(rOut, 0, zCenter + dZ).theta()
			);
  }
  return {minTheta, maxTheta};
}
*/

std::array<double, 2> envelopeEtaExtremes (uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto detelement = volMgr.lookupDetElement(aVolumeId);
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  // calculate values of eta in all possible corners of the envelope
  auto dim = envelopeDimensions(aVolumeId);
  double minEta = 0;
  double maxEta = 0;
  for (uint i = 0; i < 8; i++) {
    // coefficients to get all combinations of corners
    int iX = -1 + 2 * ((i / 2) % 2);
    int iY = -1 + 2 * (i % 2);
    int iZ = -1 + 2 * (i / 4);
    double outDimGlobal[3];
    double inDimLocal[] = {iX * dim.x(), iY * dim.y(), iZ * dim.z()};
    transformMatrix.LocalToMaster(inDimLocal, outDimGlobal);
    double eta = CLHEP::Hep3Vector(outDimGlobal[0], outDimGlobal[1], outDimGlobal[2]).eta();
    if (i == 0) {
      minEta = eta;
      maxEta = eta;
    }
    if (eta < minEta) {
      minEta = eta;
    }
    if (eta > maxEta) {
      maxEta = eta;
    }
  }
  return {minEta, maxEta};
}

// GM: unnneded, kept for debug, will remove before merge
/*
std::array<double, 2> envelopeThetaExtremes (uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto detelement = volMgr.lookupDetElement(aVolumeId);
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  // calculate values of theta in all possible corners of the envelope
  auto dim = envelopeDimensions(aVolumeId);
  double minTheta = 0;
  double maxTheta = 0;
  for (uint i = 0; i < 8; i++) {
    // coefficients to get all combinations of corners
    int iX = -1 + 2 * ((i / 2) % 2);
    int iY = -1 + 2 * (i % 2);
    int iZ = -1 + 2 * (i / 4);
    double outDimGlobal[3];
    double inDimLocal[] = {iX * dim.x(), iY * dim.y(), iZ * dim.z()};
    transformMatrix.LocalToMaster(inDimLocal, outDimGlobal);
    double theta = CLHEP::Hep3Vector(outDimGlobal[0], outDimGlobal[1], outDimGlobal[2]).theta();
    if (i == 0) {
      minTheta = theta;
      maxTheta = theta;
    }
    if (theta < minTheta) {
      minTheta = theta;
    }
    if (theta > maxTheta) {
      maxTheta = theta;
    }
  }
  // std::cout << "DEBUG: minTheta = " << minTheta << " maxTheta = " << maxTheta << std::endl;
  return {minTheta, maxTheta};
}
*/

std::array<double, 2> volumeEtaExtremes(uint64_t aVolumeId) {
  // try if volume is a cylinder/disc
  auto etaExtremes = tubeEtaExtremes(aVolumeId);
  if (etaExtremes[0] != 0 or etaExtremes[1] != 0) {
    return etaExtremes;
  } else {
    return envelopeEtaExtremes(aVolumeId);
  }
}

// GM: unnneded, kept for debug, will remove before merge
/*
std::array<double, 2> volumeThetaExtremes(uint64_t aVolumeId) {
  // try if volume is a cylinder/disc
  auto thetaExtremes = tubeThetaExtremes(aVolumeId);
  if (thetaExtremes[0] != 0 or thetaExtremes[1] != 0) {
    return thetaExtremes;
  } else {
    return envelopeThetaExtremes(aVolumeId);
  }
}
*/

std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridPhiEta& aSeg) {
  // get segmentation number of bins in phi
  uint phiCellNumber = aSeg.phiBins();
  // get segmentation cell width in eta
  double etaCellSize = aSeg.gridSizeEta();
  // get min and max eta of the volume
  auto etaExtremes = volumeEtaExtremes(aVolumeId);
  // calculate the number of eta volumes
  // max - min = full eta range, - size = not counting the middle cell centred at 0, + 1 to account for that cell
  uint cellsEta = ceil(( etaExtremes[1] - etaExtremes[0] - etaCellSize ) / 2 / etaCellSize) * 2 + 1;
  uint minEtaID = int(floor((etaExtremes[0] + 0.5 * etaCellSize - aSeg.offsetEta()) / etaCellSize));
  return {phiCellNumber, cellsEta, minEtaID};
}

std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridPhiTheta& aSeg) {
  uint phiCellNumber = aSeg.phiBins();
  double thetaCellSize = aSeg.gridSizeTheta();
  auto etaExtremes = volumeEtaExtremes(aVolumeId);
  double thetaMin = 2.*atan(exp(-etaExtremes[1]));
  double thetaMax = 2.*atan(exp(-etaExtremes[0]));

  // GM debug - remove before merge
  /*
  auto thetaExtremes = volumeThetaExtremes(aVolumeId);
  std::cout << thetaMin << " " << thetaExtremes[0] << std::endl;
  std::cout << thetaMax << " " << thetaExtremes[1] << std::endl;
  assert(abs(thetaMin-thetaExtremes[0])<1e-7);
  assert(abs(thetaMax-thetaExtremes[1])<1e-7);
  */
  uint cellsTheta = ceil(( thetaMax - thetaMin - thetaCellSize ) / 2 / thetaCellSize) * 2 + 1;
  uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - aSeg.offsetTheta()) / thetaCellSize));
  return {phiCellNumber, cellsTheta, minThetaID};
}

std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged& aSeg) {

  const dd4hep::DDSegmentation::BitFieldCoder* aDecoder = aSeg.decoder();
  int nLayer = aDecoder->get(aVolumeId, aSeg.fieldNameLayer());
  // + 0.5 to avoid integer division
  uint moduleCellNumber = ceil((aSeg.nModules() + 0.5) / aSeg.mergedModules(nLayer));
  // get minimum and maximum theta of volume
  auto etaExtremes = volumeEtaExtremes(aVolumeId);
  double thetaMin = 2.*atan(exp(-etaExtremes[1]));
  double thetaMax = 2.*atan(exp(-etaExtremes[0]));
  // GM debug - remove before merge
  /*
  auto thetaExtremes = volumeThetaExtremes(aVolumeId);
  std::cout << thetaMin << " " << thetaExtremes[0] << std::endl;
  std::cout << thetaMax << " " << thetaExtremes[1] << std::endl;
  assert(abs(thetaMin - thetaExtremes[0]) < 1e-7);
  assert(abs(thetaMax - thetaExtremes[1]) < 1e-7);
  */

  // convert to minimum and maximum theta bins
  double thetaCellSize = aSeg.gridSizeTheta();
  double thetaOffset = aSeg.offsetTheta();
  uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
  uint maxThetaID = int(floor((thetaMax + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
  // correct minThetaID and maxThetaID for merging
  uint mergedThetaCells = aSeg.mergedThetaCells(nLayer);
  minThetaID -= (minThetaID % mergedThetaCells);
  maxThetaID -= (maxThetaID % mergedThetaCells);
  uint cellsTheta = 1 + (maxThetaID - minThetaID)/ mergedThetaCells;
  // GM debug - remove before merge
  //std::cout << "DEBUG: layer = " << nLayer << " moduleCellNumber = " << moduleCellNumber << " cellsTheta = " << cellsTheta << " minThetaID = " << minThetaID << std::endl;
  return {moduleCellNumber, cellsTheta, minThetaID};
}

std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::PolarGridRPhi& aSeg) {
  // get half-widths,
  auto tubeSizes = tubeDimensions(aVolumeId);
  // get segmentation cell width
  double rCellSize = aSeg.gridSizeR();
  double phiCellSize = aSeg.gridSizePhi();
  uint cellsRout = ceil(tubeSizes.y() / rCellSize);
  uint cellsRin = floor(tubeSizes.x() / rCellSize);
  uint cellsR = cellsRout - cellsRin;
  uint cellsPhi = ceil(2 * M_PI / phiCellSize);
  return {cellsR, cellsPhi};
}

unsigned int countPlacedVolumes(TGeoVolume* aHighestVolume, const std::string& aMatchName) {
  int numberOfPlacedVolumes = 0;
  TGeoNode* node;
  TGeoIterator next(aHighestVolume);
  while ((node = next())) {
    std::string currentNodeName = node->GetName();
    if (currentNodeName.find(aMatchName) != std::string::npos) {
      ++numberOfPlacedVolumes;
    }
  }
  return numberOfPlacedVolumes;
}
}
}
