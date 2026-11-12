// TrackingData.h — value type passed through EventContext.
// Separate from Tracking.h so CollectionWriter.h can include this
// without pulling in the writer class (avoids include cycle since
// TrackingWriter inherits from CollectionWriter).

#pragma once

#include <edm4hep/MutableReconstructedParticle.h>

#include <cstdint>
#include <vector>

namespace delphi_edm4hep::tracking {

// How a <prefix>_MAIN_Particles entry relates to SKELANA's VECP list,
// published per particle as <prefix>_MAIN_Particles_VecpStatus.
//
// SKELANA (IFLODR=1, IFLRVR=111) removes the charged VECP row of a track that
// fails the charged selection and is not rescued by the beam-spot refit or
// the Mammoth step, but whose calorimeter energy passes the neutral cuts
// (PSHCTRECCAL). The same PA is then re-entered as a NEUTRAL VECP row by
// PSHCHREC: LVECP = the original LPA, VECP(7) = 0, a massless 4-vector of
// the passing HPC+EMF+HCAL+STIC energy along the first good shower.
// The PA keeps its own entry (it owns the track); the recovered 4-vector
// becomes a second, neutral entry.
enum class VecpStatus : std::int32_t {
  OwnRow           = 0,  // 4-vector and verdict from the PA's own VECP row
  RecoveredSource  = 1,  // charged PA recovered as neutral by SKELANA:
                         //   p, E, m are NaN; see the RecoveredNeutral entry
  RecoveredNeutral = 2,  // neutral 4-vector SKELANA built from a charged
                         //   PA's calorimeter deposits; particles[0] is the
                         //   RecoveredSource entry
  NoRow            = 3,  // charged PA with no VECP row at all: p, E, m NaN
};

// Output of TrackingWriter::emit():
//   - particle_handles: one MutableReconstructedParticle handle per
//     emitted Particle, in PA-walk order (a RecoveredNeutral entry directly
//     follows its source). Handles remain valid after the collection has
//     been moved into the Frame.
//   - vecp_to_particle: VECP index (1..NVECP, with sk::LVPART=1) ->
//     0-based index into the emitted Particle collection. -1 if no
//     Particle was emitted for that VECP entry. Index 0 unused.
//   - pa_to_particle:   PA-walk index (0-based) -> 0-based Particle
//     collection index, or -1 if Tracking dropped the PA. For a recovered
//     PA this is the track-owning RecoveredSource entry.
//   - pa_to_recovered_neutral: PA-walk index -> 0-based index of the
//     RecoveredNeutral entry built from that PA, or -1. The PA's showers
//     belong to this entry, since its energy is made from them.
struct Output {
  std::vector<edm4hep::MutableReconstructedParticle> particle_handles;
  std::vector<int> vecp_to_particle;
  std::vector<int> pa_to_particle;
  std::vector<int> pa_to_recovered_neutral;

  // The Particle that owns this PA's calorimeter showers, or -1.
  int calorimeterOwner(int paIdx) const {
    if (paIdx < 0) return -1;
    if (paIdx < static_cast<int>(pa_to_recovered_neutral.size()) &&
        pa_to_recovered_neutral[paIdx] >= 0) {
      return pa_to_recovered_neutral[paIdx];
    }
    return paIdx < static_cast<int>(pa_to_particle.size())
               ? pa_to_particle[paIdx] : -1;
  }
};

}  // namespace delphi_edm4hep::tracking
