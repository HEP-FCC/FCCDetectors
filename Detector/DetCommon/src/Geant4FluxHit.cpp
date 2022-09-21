
#include "DetCommon/Geant4FluxHit.h"

namespace k4 {

// G4 allocation method
G4ThreadLocal G4Allocator<Geant4FluxHit>* Geant4FluxHitAllocator = 0;
// Destructor
Geant4FluxHit::~Geant4FluxHit() {}
// Default Constructor
Geant4FluxHit::Geant4FluxHit() {}
// Constructor setting some members
Geant4FluxHit::Geant4FluxHit(unsigned int aTrackId,
                             int aPdgId,
                             double aParticleFlux,
                             double aEnergy,
                             double aTime) : trackId(aTrackId),
                                             pdgId(aPdgId),
                                             particleFlux(aParticleFlux),
                                             energy(aEnergy),
                                             time(aTime) {}

// comparison operator
G4int Geant4FluxHit::operator==(const Geant4FluxHit& right) const {
  return (this == &right) ? 1 : 0;
}

}  // namespace k4
