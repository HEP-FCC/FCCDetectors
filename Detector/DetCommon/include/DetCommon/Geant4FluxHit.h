#ifndef SIMG4COMMON_GEANT4FLUXHIT
#define SIMG4COMMON_GEANT4FLUXHIT

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"

// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace k4 {

/** @class  Geant4FluxHit
 *
 * Data structure to hold the geant4 output of the Fluxmeter.
 *
 */
class Geant4FluxHit : public G4VHit {
public:
  /// Default constructor
  Geant4FluxHit();
  /// Constructor setting some members
  Geant4FluxHit(unsigned int aTrackId,
                int aPdgId,
                double aParticleFlux,
                double aEnergy,
                double aTime);
  // Destructor
  virtual ~Geant4FluxHit();

  /// Comparison operator
  G4int operator==(const Geant4FluxHit&) const;
  /// New operator needed for g4 memory allocation
  inline void* operator new(size_t);
  /// Delete operator needed for g4 memory allocation
  inline void operator delete(void*);

  /// Method from base class, unused
  virtual void Draw(){};
  /// Method from base class, unused
  virtual void Print(){};

  // These members are public, following the example of G4VHit:

  /// The pre-step position of the step in which energy was deposited
  CLHEP::Hep3Vector position;
  /// The DD4hep cellID of the volume in which the energy was deposited
  unsigned long cellId;
  /// The g4 trackId of the particle that deposited the energy
  unsigned int trackId;
  /// The particle data group identification code for the particle
  int pdgId;
  /// The particle flux into the fluxmeter
  double particleFlux;
  /// The particle kinetic energy into the fluxmeter
  double energy;
  /// The particle total energy into the fluxmeter
  double totalEnergy;
  /// The time coordinate of the energy deposit
  double time;
  /// The creation vertex of the incident particle
  CLHEP::Hep3Vector particleVertex;
  /// The momentum of the incident particle when entering the fluxmeter
  CLHEP::Hep3Vector momentum;
  /// The particle charge when entering the fluxmeter
  double charge;
  /// The particle mass when entering the fluxmeter
  double mass;
};

// types and functions for G4 memory allocation, inspired by the G4VHit classes
// in Geant4 examples

extern G4ThreadLocal G4Allocator<Geant4FluxHit>* Geant4FluxHitAllocator;

inline void* Geant4FluxHit::operator new(size_t) {
  if (!Geant4FluxHitAllocator) {
    Geant4FluxHitAllocator = new G4Allocator<Geant4FluxHit>;
  }

  return (void*)Geant4FluxHitAllocator->MallocSingle();
}

inline void Geant4FluxHit::operator delete(void* hit) {
  Geant4FluxHitAllocator->FreeSingle((Geant4FluxHit*)hit);
}

}  // namespace k4

#endif /* SIMG4COMMON_GEANT4FLUXHIT */
