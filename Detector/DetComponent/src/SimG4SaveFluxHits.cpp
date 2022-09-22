#include "SimG4SaveFluxHits.h"

#include "DetCommon/Geant4FluxHit.h"
#include "DetCommon/Units.h"

// k4SimGeant4
#include "k4Interface/IGeoSvc.h"

// Geant4
#include "G4Event.hh"

// datamodel
#include "edm4hep/Vector3d.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Segmentations.h"


DECLARE_COMPONENT(SimG4SaveFluxHits)

SimG4SaveFluxHits::SimG4SaveFluxHits(const std::string& aType, const std::string& aName, const IInterface* aParent)
    : GaudiTool(aType, aName, aParent), m_geoSvc("GeoSvc", aName), m_eventDataSvc("EventDataSvc", "SimG4SaveCalHits")  {
  declareInterface<ISimG4SaveOutputTool>(this);
  declareProperty("HitPosition", m_position, "Handle for hit position");
  declareProperty("HitCellID", m_cellID, "Handle hit cell ID");
  declareProperty("TrackID", m_trackID, "Handle for track ID of the incident particle");
  declareProperty("ParticlePDGID", m_pdgID, "Handle for PDG ID of the incident particle");
  declareProperty("ParticleFlux", m_particleFlux, "Handle for particle flux");
  declareProperty("ParticleEnergy", m_energy, "Handle for energy of incoming particle");
  declareProperty("HitTime", m_time, "Handle for hit time");
  declareProperty("ParticleVertex", m_particleVertex, "Handle for creation vertex of the incident particle");
  declareProperty("ParticleMomentum", m_momentum, "Handle for momentum of the particle when entering fluxmeter");
  declareProperty("ParticleCharge", m_charge, "Handle for charge at hit location");
  declareProperty("ParticleMass", m_mass, "Handle for mass at the hit location");
  declareProperty("GeoSvc", m_geoSvc);
}

SimG4SaveFluxHits::~SimG4SaveFluxHits() {}

StatusCode SimG4SaveFluxHits::initialize() {
  if (GaudiTool::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  auto lcdd = m_geoSvc->lcdd();
  auto allReadouts = lcdd->readouts();
  for (auto& readoutName : m_readoutNames) {
    if (allReadouts.find(readoutName) == allReadouts.end()) {
      error() << "Readout " << readoutName << " not found! Please check tool configuration." << endmsg;
      return StatusCode::FAILURE;
    } else {
      debug() << "Hits will be saved to EDM from the collection " << readoutName << endmsg;
    }
  }

  StatusCode sc = m_eventDataSvc.retrieve();
  m_podioDataSvc = dynamic_cast<PodioDataSvc*>(m_eventDataSvc.get());
  if (sc == StatusCode::FAILURE) {
    error() << "Error retrieving Event Data Service" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode SimG4SaveFluxHits::finalize() { return GaudiTool::finalize(); }

StatusCode SimG4SaveFluxHits::saveOutput(const G4Event& aEvent) {
  G4HCofThisEvent* collections = aEvent.GetHCofThisEvent();
  G4VHitsCollection* collect;
  k4::Geant4FluxHit* hit;
  if (collections != nullptr) {
    auto positionVec = m_position.createAndPut();
    auto cellIDVec = m_cellID.createAndPut();
    auto trackIDVec = m_trackID.createAndPut();
    auto pdgIDVec = m_pdgID.createAndPut();
    auto particleFluxVec = m_particleFlux.createAndPut();
    auto energyVec = m_energy.createAndPut();
    auto timeVec = m_time.createAndPut();
    auto particleVertexVec = m_particleVertex.createAndPut();
    auto momentumVec = m_momentum.createAndPut();
    auto chargeVec = m_charge.createAndPut();
    auto massVec = m_mass.createAndPut();

    for (int iter_coll = 0; iter_coll < collections->GetNumberOfCollections(); iter_coll++) {
      collect = collections->GetHC(iter_coll);
      if (std::find(m_readoutNames.begin(), m_readoutNames.end(), collect->GetName()) != m_readoutNames.end()) {
        // Add CellID encoding string to collection metadata
        // auto lcdd = m_geoSvc->lcdd();
        // auto allReadouts = lcdd->readouts();
        // auto idspec = lcdd->idSpecification(collect->GetName());
        // auto field_str = idspec.fieldDescription();
        // auto& coll_md = m_podioDataSvc->getProvider().getCollectionMetaData( m_cellID.get()->getID() );
        // coll_md.setValue("CellIDEncodingString", field_str);

        size_t n_hit = collect->GetSize();
        debug() << "\t" << n_hit << " hits are stored in a collection #" << iter_coll << ": " << collect->GetName()
                << endmsg;
        for (size_t iter_hit = 0; iter_hit < n_hit; iter_hit++) {
          hit = dynamic_cast<k4::Geant4FluxHit*>(collect->GetHit(iter_hit));

          positionVec->emplace_back(std::vector<double>{
            hit->position.x() * sim::g42edm::length,
            hit->position.y() * sim::g42edm::length,
            hit->position.z() * sim::g42edm::length,
          });
          cellIDVec->emplace_back(hit->cellId);
          trackIDVec->emplace_back(hit->trackId);
          pdgIDVec->emplace_back(hit->pdgId);
          particleFluxVec->emplace_back(hit->particleFlux);
          energyVec->emplace_back(hit->energy * sim::g42edm::energy);
          timeVec->emplace_back(hit->time);
          particleVertexVec->emplace_back(std::vector<double>{
            hit->particleVertex.x() * sim::g42edm::length,
            hit->particleVertex.y() * sim::g42edm::length,
            hit->particleVertex.z() * sim::g42edm::length,
          });
          momentumVec->emplace_back(std::vector<double>{
            hit->momentum.x() * sim::g42edm::energy,
            hit->momentum.y() * sim::g42edm::energy,
            hit->momentum.z() * sim::g42edm::energy,
          });
          chargeVec->emplace_back(hit->charge);
          massVec->emplace_back(hit->mass);
        }
      }
    }
  }
  return StatusCode::SUCCESS;
}
