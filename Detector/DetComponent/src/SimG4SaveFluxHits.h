#ifndef DETCOMPONENTS_G4SAVEFLUXHITS_H
#define DETCOMPONENTS_G4SAVEFLUXHITS_H

// Gaudi
#include "GaudiAlg/GaudiTool.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4Interface/ISimG4SaveOutputTool.h"
class IGeoSvc;

// datamodel
namespace edm4hep {
class Vector3d;
}

/** @class SimG4SaveFluxHits SimG4Components/src/SimG4SaveFluxHits.h SimG4SaveFluxHits.h
 *
 *  Save calorimeter hits tool.
 *  All collections passed in the job options will be saved (\b'readoutNames').
 *  Readout name is defined in DD4hep XML file as the attribute 'readout' of 'detector' tag.
 *  If (\b'readoutNames') contain no elements or names that do not correspond to any hit collection,
 *  tool will fail at initialization.
 *  [For more information please see](@ref md_sim_doc_geant4fullsim).
 *
 *  @author Anna Zaborowska
 */

class SimG4SaveFluxHits : public GaudiTool, virtual public ISimG4SaveOutputTool {
public:
  explicit SimG4SaveFluxHits(const std::string& aType,
                             const std::string& aName,
                             const IInterface* aParent);
  virtual ~SimG4SaveFluxHits();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize();
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize();
  /**  Save the data output.
   *   Saves the calorimeter hits from the collections as specified in the job options in \b'readoutNames'.
   *   @param[in] aEvent Event with data to save.
   *   @return status code
   */
  virtual StatusCode saveOutput(const G4Event& aEvent) final;

private:
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Pointer to Podio and Event Data Services
  PodioDataSvc* m_podioDataSvc;
  ServiceHandle<IDataProviderSvc> m_eventDataSvc;
  /// The pre-step position of the step in which energy was deposited
  DataHandle<std::vector<std::vector<double>>> m_position{"HitPosition", Gaudi::DataHandle::Writer, this};
  /// The DD4hep cellID of the volume in which the energy was deposited
  DataHandle<std::vector<unsigned long>> m_cellID{"HitCellID", Gaudi::DataHandle::Writer, this};
  /// The g4 trackId of the particle that deposited the energy
  DataHandle<std::vector<unsigned int>> m_trackID{"TrackID", Gaudi::DataHandle::Writer, this};
  /// The particle data group identification code for the particle
  DataHandle<std::vector<int>> m_pdgID{"ParticlePDGID", Gaudi::DataHandle::Writer, this};
  /// The particle flux into the fluxmeter
  DataHandle<std::vector<double>> m_particleFlux{"ParticleFlux", Gaudi::DataHandle::Writer, this};
  /// The kinetic energy of incoming particle into the fluxmeter
  DataHandle<std::vector<double>> m_energy{"ParticleEnergy", Gaudi::DataHandle::Writer, this};
  /// The total energy of incoming particle into the fluxmeter
  DataHandle<std::vector<double>> m_totalEnergy{"ParticleTotalEnergy", Gaudi::DataHandle::Writer, this};
  /// The time coordinate of the hit
  DataHandle<std::vector<double>> m_time{"HitTime", Gaudi::DataHandle::Writer, this};
  /// The creation vertex of the incoming particle
  DataHandle<std::vector<std::vector<double>>> m_particleVertex{"ParticleVertex", Gaudi::DataHandle::Writer, this};
  /// The momentum of the incoming particle when entering fluxmeter
  DataHandle<std::vector<std::vector<double>>> m_momentum{"ParticleMomentum", Gaudi::DataHandle::Writer, this};
  /// The particle charge when entering the fluxmeter
  DataHandle<std::vector<double>> m_charge{"ParticleCharge", Gaudi::DataHandle::Writer, this};
  /// The particle mass when entering the fluxmeter
  DataHandle<std::vector<double>> m_mass{"ParticleMass", Gaudi::DataHandle::Writer, this};
  /// Name of the readouts (hits collections) to save
  Gaudi::Property<std::vector<std::string>> m_readoutNames{
      this, "readoutNames", {}, "Name of the readouts (hits collections) to save"};
};

#endif /* DETCOMPONENTS_G4SAVEFLUXHITS_H */
