#include "DetSensitive/FluxMeter.h"

// FCCDetectors
#include "DetCommon/DetUtils.h"
#include "DetCommon/Geant4FluxHit.h"

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
  FluxMeter::FluxMeter(const std::string& aDetectorName,
                       const std::string& aReadoutName,
                       const dd4hep::Segmentation& aSeg)
      : G4VSensitiveDetector(aDetectorName),
        m_fluxCollection(nullptr),
        m_seg(aSeg) {
    // name of the collection of hits is determined by the readout name (from XML)
    collectionName.insert(aReadoutName);
  }

  FluxMeter::~FluxMeter() {}

  void FluxMeter::Initialize(G4HCofThisEvent* aHitsCollections) {
    // create a collection of hits and add it to G4HCofThisEvent
    // deleted in ~G4Event
    m_fluxCollection = new G4THitsCollection<k4::Geant4FluxHit>(SensitiveDetectorName,
                                                                collectionName[0]);
    aHitsCollections->AddHitsCollection(
      G4SDManager::GetSDMpointer()->GetCollectionID(m_fluxCollection),
      m_fluxCollection
    );

    // Determining if detector is a tube or box
    // auto fluxDet = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName);
    //fluxDet

  }

  bool FluxMeter::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
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
    uint64_t cellID = utils::cellID(m_seg, *aStep);
    // std::vector<double> cellDim = m_seg.segmentation()->cellDimensions(cellID);
    double cellArea = 1.;

    // Getting track info
    G4Track* aTrack = aStep->GetTrack();

    // Flux into the flux meter
    double particleFlux = 1 * preStep->GetWeight()
                            * angleFactor
                            / cellArea;

    auto hit = new k4::Geant4FluxHit(
      aTrack->GetTrackID(),
      aTrack->GetDefinition()->GetPDGEncoding(),
      particleFlux,
      preStep->GetTotalEnergy(),
      aTrack->GetGlobalTime()
    );
    hit->cellId = cellID;
    hit->position = preStep->GetPosition();
    hit->particleVertex = aTrack->GetVertexPosition();
    hit->momentum = preStep->GetMomentum();
    hit->mass = preStep->GetMass();
    hit->charge = preStep->GetCharge();
    m_fluxCollection->insert(hit);

    // Should the particle be killed or left alone?
    // aTrack->SetTrackStatus(fStopAndKill);
    // return true;

    return true;
  }
}
