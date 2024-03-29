#ifndef DETSENSITIVE_AGGREGATECALORIMETERSD_H
#define DETSENSITIVE_AGGREGATECALORIMETERSD_H


// DD4hep
#include "DD4hep/Segmentations.h"

// FCCSW
#include "DetCommon/Geant4CaloHit.h"

// Geant
#include "G4THitsCollection.hh"
#include "G4VSensitiveDetector.hh"

namespace k4 {
class Geant4CaloHit;
}

/** AggregateCalorimeterSD DetectorDescription/DetSensitive/src/AggregateCalorimeterSD.h AggregateCalorimeterSD.h
 *
 *  Sensitive detector for calorimeter (aggregates energy deposits within each cell).
 *  It is based on dd4hep::sim::Geant4GenericSD<Calorimeter> (but it is not identical).
 *  In particular, the position of the hit is set to G4Step::GetPreStepPoint() position.
 *  No timing information is saved (energy deposits are aggregated in the cells)
 *
 *  @author    Anna Zaborowska
 */

namespace det {
class AggregateCalorimeterSD : public G4VSensitiveDetector {
public:
  /** Constructor.
   *  @param aDetectorName Name of the detector
   *  @param aReadoutName Name of the readout (used to name the collection)
   *  @param aSeg Segmentation of the detector (used to retrieve the cell ID)
   */
  AggregateCalorimeterSD(const std::string& aDetectorName,
                         const std::string& aReadoutName,
                         const dd4hep::Segmentation& aSeg);
  /// Destructor
  virtual ~AggregateCalorimeterSD();
  /** Initialization.
   *  Creates the hit collection with the name passed in the constructor.
   *  The hit collection is registered in Geant.
   *  @param aHitsCollections Geant hits collection.
   */
  virtual void Initialize(G4HCofThisEvent* aHitsCollections) final;
  /** Process hit once the particle hit the sensitive volume.
   *  Checks if the energy deposit is larger than 0, calculates the position and cellID,
   *  saves that into the hit collection.
   *  If there is already entry in the same cell, the energy is accumulated.
   *  Otherwise new hit is created.
   *  @param aStep Step in which particle deposited the energy.
   */
  virtual bool ProcessHits(G4Step* aStep, G4TouchableHistory*) final;

private:
  /// Collection of calorimeter hits
  G4THitsCollection<k4::Geant4CaloHit>* m_calorimeterCollection;
  /// Segmentation of the detector used to retrieve the cell Ids
  dd4hep::Segmentation m_seg;
};
}

#endif /* DETSENSITIVE_AGGREGATECALORIMETERSD_H */
