#include "DetSensitive/FluxSD.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetCommon/Geant4CaloHit.h"

// DD4hep
#include "DDG4/Defs.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4VolumeManager.h"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// Geant4
#include "G4SDManager.hh"
#include "G4VSolid.hh"


namespace det {
  FluxSD::FluxSD(const std::string& aDetectorName,
                 const std::string& aReadoutName,
                 const dd4hep::Segmentation& aSeg)
      : G4VSensitiveDetector(aDetectorName),
        m_fluxCollection(nullptr),
        m_seg(aSeg) {
    // name of the collection of hits is determined by the readout name (from XML)
    collectionName.insert(aReadoutName);
  }

  FluxSD::~FluxSD() {}

  void FluxSD::Initialize(G4HCofThisEvent* aHitsCollections) {
    // create a collection of hits and add it to G4HCofThisEvent
    // deleted in ~G4Event
    m_fluxCollection = new G4THitsCollection<k4::Geant4CaloHit>(SensitiveDetectorName,
                                                                collectionName[0]);
    aHitsCollections->AddHitsCollection(
      G4SDManager::GetSDMpointer()->GetCollectionID(m_fluxCollection),
      m_fluxCollection
    );

    // Determining if detector is a tube or box
    // auto fluxDet = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName);
    //fluxDet

  }

  bool FluxSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4StepPoint* preStep = aStep->GetPreStepPoint();
    // Determine if step is on the boundary
    // I'm not sure if check on boundary closeness is needed
    if (preStep->GetStepStatus() != fGeomBoundary) {
      return true;
    }

    // Get local position of the prestep
    G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
    G4ThreeVector global = preStep->GetPosition();
    G4ThreeVector local =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(global);

    // Get surface normal to the surface
    G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
    G4VSolid* solid = physVol->GetLogicalVolume()->GetSolid();
    G4ThreeVector surfNorm = solid->SurfaceNormal(local);

    // Get local direction of the passing particle
    G4ThreeVector globDir = preStep->GetMomentumDirection();
    G4ThreeVector localDir =
      theTouchable->GetHistory()->GetTopTransform().TransformAxis(globDir);

    // Calculate angle factor
    double angleFactor = surfNorm.dot(localDir);
    if (angleFactor < 0) {
      angleFactor *= -1;
    }
    // std::cout << "Angle factor: " << angleFactor << std::endl;

    // Get cell area (Not used, needs correct segmentation)
    uint64_t cellId = utils::cellID(m_seg, *aStep);
    // std::vector<double> cellDim = m_seg.segmentation()->cellDimensions(cellId);
    double cellArea = 1.;

    // Flux into the flux meter
    double flux = 1 * preStep->GetWeight()
                    * angleFactor
                    / cellArea;
    // std::cout << "flux: " << flux << std::endl;

    // Getting track info
    // G4Track* aTrack = aStep->GetTrack();

    // Should the particle be killed or left alone?
    // aTrack->SetTrackStatus(fStopAndKill);
    // return true;

    // Part from dd4hep::sim::Geant4GenericSD<Calorimeter>
    CLHEP::Hep3Vector prePos = aStep->GetPreStepPoint()->GetPosition();
    CLHEP::Hep3Vector postPos = aStep->GetPostStepPoint()->GetPosition();
    CLHEP::Hep3Vector midPos = 0.5 * (postPos + prePos);

    k4::Geant4CaloHit* hit = nullptr;
    k4::Geant4CaloHit* hitMatch = nullptr;
    // Check if there is already some energy deposit in that cell
    for (size_t i = 0; i < m_fluxCollection->entries(); i++) {
      hit = dynamic_cast<k4::Geant4CaloHit*>(m_fluxCollection->GetHit(i));
      if (hit->cellID == cellId) {
        hitMatch = hit;
        hitMatch->energyDeposit += flux;
        return true;
      }
    }
    // if not, create a new hit
    // deleted in ~G4Event
    hitMatch = new k4::Geant4CaloHit(0, // track->GetTrackID()
                                     0, // track->GetDefinition()->GetPDGEncoding()
                                     flux,
                                     0 // track ->GetGlobalTime()
                                     ) ;

    hitMatch->position = midPos;
    hitMatch->cellID = cellId;
    m_fluxCollection->insert(hitMatch);

    return true;
  }
}
