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
						   bool aDiagonal) {
  std::vector<uint64_t> neighbours;
  int nLayer = aDecoder.get(aCellId, aSeg.fieldNameLayer());
  dd4hep::DDSegmentation::CellID cID = aCellId;

  // find index of module in field extremes vector
  int idModuleField = -1;
  for (uint itField = 0; itField < aFieldNames.size(); itField++) {
    if (aFieldNames[itField] == aSeg.fieldNameModule()) {
      idModuleField = (int) itField;
      break;
    }
  }
  if (idModuleField < 0) {
    std::cout << "WARNING: module field " << aSeg.fieldNameModule() << " not found in aFieldNames vector" << std::endl;
    std::cout << "WARNING: will return empty neighbour map" << std::endl;
    return neighbours;
  }

  // now find the neighbours
  for (uint itField = 0; itField < aFieldNames.size(); itField++) {
    // retrieve name of field and corresponding value in cellID
    const auto& field = aFieldNames[itField];
    //dd4hep::DDSegmentation::CellID id = aDecoder.get(aCellId, field);
    int id = aDecoder.get(aCellId, field);

    if (field == aSeg.fieldNameLayer()) {

      // for neighbours across different layers, we have to take into
      // account that merging along module and/or theta could be different
      // so one cell in layer N could be neighbour to several in layer N+-1
      int module_id = aDecoder.get(aCellId, aSeg.fieldNameModule());
      int theta_id = aDecoder.get(aCellId, aSeg.fieldNameTheta());
      
      // The cells are classified in a different way whether they are
      // direct neighbours (common surface), diagonal neighbours (common edge or vertex)
      // or neither.
      // To decide this, we need to check how the cells are related in both directions:
      // neighbours (edge at least partially in common), diagonal neigbours (common vertex),
      // none
      int neighbourTypeModuleDir; // 0: not neighbour; 1: diagonal neighbour; 2: neighbour in module direction
      int neighbourTypeThetaDir; // 0: not neighbour; 1: diagonal neighbour; 2: neighbour in module direction

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

	// OLD VERSION
	/*
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
	*/

	// NEW VERSION
	for (int i=-1; i <= aSeg.mergedThetaCells(nLayer); i++) {
	  int theta_id_neighbour = (theta_id + i) - ((theta_id + i) % aSeg.mergedThetaCells(nLayer+deltaLayer));
	  if (theta_id_neighbour >= theta_id && theta_id_neighbour < (theta_id + aSeg.mergedThetaCells(nLayer))) neighbourTypeThetaDir = 2;
	  else if (theta_id_neighbour < theta_id && theta_id_neighbour > (theta_id - aSeg.mergedThetaCells(nLayer+deltaLayer))) neighbourTypeThetaDir = 2;
	  else if (theta_id_neighbour == (theta_id + aSeg.mergedThetaCells(nLayer))) neighbourTypeThetaDir = 1;
	  else if (theta_id_neighbour == (theta_id - aSeg.mergedThetaCells(nLayer+deltaLayer))) neighbourTypeThetaDir = 1;
	  else neighbourTypeThetaDir = 0;
		  
	  // if there is no point of contact along theta, no need to check also for module direction
	  if (neighbourTypeThetaDir == 0) continue;
	  // if we are not considering diagonal neighbours, and cells in theta have only an edge in common, then skip
	  if (!aDiagonal && neighbourTypeThetaDir == 1) continue; 
	  // otherwise, check status along module direction
	  for (int j=-1; j <= aSeg.mergedModules(nLayer); j++) {
	    int module_id_neighbour = (module_id + j) - ((module_id + j) % aSeg.mergedModules(nLayer+deltaLayer));
	    int module_id_neighbour_cyclic = cyclicNeighbour(module_id_neighbour,
							     aFieldExtremes[idModuleField]
							     );

	    if (module_id_neighbour >= module_id && module_id_neighbour < (module_id + aSeg.mergedModules(nLayer))) neighbourTypeModuleDir = 2;
	    else if (module_id_neighbour < module_id && module_id_neighbour > (module_id - aSeg.mergedModules(nLayer+deltaLayer))) neighbourTypeModuleDir = 2;
	    else if (module_id_neighbour == (module_id + aSeg.mergedModules(nLayer))) neighbourTypeModuleDir = 1;
	    else if (module_id_neighbour == (module_id - aSeg.mergedModules(nLayer+deltaLayer))) neighbourTypeModuleDir = 1;
	    else neighbourTypeModuleDir = 0;
	  
	    // if there is no point of contact along module, then skip
	    if (neighbourTypeModuleDir == 0) continue;
	    // otherwise: if neighbours along both module and theta, or along one of the two
	    // and we also consider diagonal neighbours, then add cells to list of neighbours
	    if ( (neighbourTypeModuleDir == 2 && neighbourTypeThetaDir==2) || 
		 (aDiagonal && ((neighbourTypeThetaDir > 1) || (neighbourTypeModuleDir>1))) ) {
	      aDecoder.set(cID, aSeg.fieldNameModule(), module_id_neighbour_cyclic);
	      aDecoder.set(cID, aSeg.fieldNameTheta(), theta_id_neighbour);
	      neighbours.emplace_back(cID);
	    }
	  }
	}
	// end NEW VERSION
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
  /*
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
  */
  
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

std::array<double, 2> volumeEtaExtremes(uint64_t aVolumeId) {
  // try if volume is a cylinder/disc
  auto etaExtremes = tubeEtaExtremes(aVolumeId);
  if (etaExtremes[0] != 0 or etaExtremes[1] != 0) {
    return etaExtremes;
  } else {
    return envelopeEtaExtremes(aVolumeId);
  }
}

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

  uint cellsTheta = ceil(( thetaMax - thetaMin - thetaCellSize ) / 2 / thetaCellSize) * 2 + 1;
  uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - aSeg.offsetTheta()) / thetaCellSize));
  return {phiCellNumber, cellsTheta, minThetaID};
}

std::array<uint, 3> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged& aSeg) {

  const dd4hep::DDSegmentation::BitFieldCoder* aDecoder = aSeg.decoder();
  int nLayer = aDecoder->get(aVolumeId, aSeg.fieldNameLayer());
  // + 0.5 to avoid integer division
  uint nModules = aSeg.nModules() / aSeg.mergedModules(nLayer);
  if (aSeg.nModules() % aSeg.mergedModules(nLayer) != 0) nModules++;
  // get minimum and maximum theta of volume
  auto etaExtremes = volumeEtaExtremes(aVolumeId);
  double thetaMin = 2.*atan(exp(-etaExtremes[1]));
  double thetaMax = 2.*atan(exp(-etaExtremes[0]));

  // convert to minimum and maximum theta bins
  double thetaCellSize = aSeg.gridSizeTheta();
  double thetaOffset = aSeg.offsetTheta();
  uint minThetaID = int(floor((thetaMin + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
  uint maxThetaID = int(floor((thetaMax + 0.5 * thetaCellSize - thetaOffset) / thetaCellSize));
  // correct minThetaID and maxThetaID for merging
  uint mergedThetaCells = aSeg.mergedThetaCells(nLayer);
  minThetaID -= (minThetaID % mergedThetaCells);
  maxThetaID -= (maxThetaID % mergedThetaCells);
  uint nThetaCells = 1 + (maxThetaID - minThetaID)/ mergedThetaCells;
  return {nModules, nThetaCells, minThetaID};
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
