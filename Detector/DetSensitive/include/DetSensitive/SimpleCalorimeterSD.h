#ifndef DETSENSITIVE_SIMPLECALORIMETERSD_H
#define DETSENSITIVE_SIMPLECALORIMETERSD_H

// DD4hep
#include "DD4hep/Segmentations.h"

// Geant
#include "G4THitsCollection.hh"
#include "G4VSensitiveDetector.hh"

namespace k4 {
class Geant4CaloHit;

}

/** SimpleCalorimeterSD DetectorDescription/DetSensitive/src/SimpleCalorimeterSD.h SimpleCalorimeterSD.h
 *
 *  Simple sensitive detector for calorimeter.
 *  It is based on dd4hep::sim::Geant4GenericSD<Calorimeter> (but it is not identical).
 *  In particular, the position of the hit is set to G4Step::GetPreStepPoint() position.
 *  New hit is created for each energy deposit.
 *  No timing information is saved.
 *
 *  @author    Anna Zaborowska
 */

namespace det {
class SimpleCalorimeterSD : public G4VSensitiveDetector {
public:
  /** Constructor.
   *  @param aDetectorName Name of the detector
   *  @param aReadoutName Name of the readout (used to name the collection)
   *  @param aSeg Segmentation of the detector (used to retrieve the cell ID)
   */
  SimpleCalorimeterSD(const std::string& aDetectorName,
                      const std::string& aReadoutName,
                      const dd4hep::Segmentation& aSeg);
  /// Destructor
  virtual ~SimpleCalorimeterSD();
  /** Initialization.
   *  Creates the hit collection with the name passed in the constructor.
   *  The hit collection is registered in Geant.
   *  @param aHitsCollections Geant hits collection.
   */
  virtual void Initialize(G4HCofThisEvent* aHitsCollections) final;
  /** Process hit once the particle hit the sensitive volume.
   *  Checks if the energy deposit is larger than 0, calculates the position and cellID,
   *  saves that into the hit collection.
   *  New hit is created for each energy deposit.
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

#endif /* DETSENSITIVE_SIMPLECALORIMETERSD_H */
