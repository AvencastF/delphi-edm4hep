// Tracking domain: PA.MAIN + PA.TRAC + PSCVEC (VECP / LVLOCK).
// TrackingWriter emits:
//   <prefix>_TRAC_Tracks                   Track + AtIP TrackState + 5x5 cov
//   <prefix>_MAIN_Particles                charged + neutral, VECP 4-mom;
//                                    NaN 4-mom for a charged PA without a
//                                    charged VECP row, plus SKELANA's
//                                    recovered neutral as its own entry
//   <prefix>_VECP_Particles_SelectionFlag  UserData<int32>
//   <prefix>_VECP_Particles_SelectionWord  UserData<int32>, SKELANA's LVSELE
//                                    of the entry's own row, raw; -1 without
//                                    one; a recovered neutral gets its own
//                                    recovery's work word (Tracking.cpp)
//   <prefix>_MAIN_Particles_VecpStatus     UserData<int32>, VecpStatus
//                                    (TrackingData.h)
//   <prefix>_MAIN_Particles_ChargeCode     UserData<int32>, raw PA.MAIN
//                                    charge code (0/1/2/3)
//   <prefix>_VECP_Particles_Charge         UserData<float>, SKELANA's VECP(7)
//                                    of the entry's own row; NaN without one
//   <prefix>_MAIN_Particles_IdentificationWord / _MassCode / _EnergyError
//                                    raw PA words IQ(LPA+3) (int32), PA.MAIN
//                                    word 10 (int32) and word 19 (float)
//   <prefix>_QTRAC_Tracks_Perigee*   UserData<float>, QTRAC(4..23) raw
//                                    (DELPHI units), parallel to TRAC_Tracks;
//                                    NaN without an own VECP row
//   <prefix>_QTRAC_Tracks_PerigeeRefitted
//                                    UserData<int32>: 1 = QTRAC(4..23) differ
//                                    from PA.TRAC (SKELANA refit), 0 = equal,
//                                    -1 = no own VECP row
//   <prefix>_QTRAC_Tracks_FitQuality / _EnergyError
//                                    UserData<float>, QTRAC(28) / QTRAC(41)
//                                    raw: PA.MAIN words 16 / 19 unless
//                                    refitted, then CONFPV's chi2 probability
//                                    / the refit momentum error (GeV); NaN
//                                    without an own VECP row
//   <prefix>_QTRAC_Tracks_d0PV / _z0PV / _d0BS
//                                    UserData<float> from QTRAC(38..40); mm
//                                    (cm x10), DELPHI sign, parallel to
//                                    TRAC_Tracks (charged-only), filled in
//                                    data and MC; SKELANA's -999 -> NaN.
//                                    The geometric <prefix>_PV_Tracks_d0PV
//                                    (Vertex.cpp; mm, LCIO sign) has
//                                    d0 ~= -1x this, z0 agrees.
//
// A companion array is named after the collection it parallels.
// and stores its Output (handle list + index maps) into ctx_.tracking
// for downstream writers.

#pragma once

#include "delphi_edm4hep/CollectionWriter.h"
#include "delphi_edm4hep/Tracking/TrackingData.h"  // Output

#include <edm4hep/TrackState.h>

#include <cstdint>

namespace delphi_edm4hep::tracking {

/// Bits of `MAIN_Particles_DetectorMask`, as masks. DELPHI numbers these bits
/// from 1, so a documented "bit N" is `1 << (N - 1)`. Read the C way every cut
/// lands one detector off and still returns plausible numbers.
///
/// Validated on 1994 data: over high-momentum barrel tracks #kDetTpc is set on
/// ~96% and #kDetId on ~73%, so an ID fraction above the TPC fraction means the
/// numbering is off by one. #kDetOd appears only above theta ~42 deg and
/// #kDetFcb only below ~35 deg, which is those detectors' acceptance.
/// @collection{MAIN_Particles_DetectorMask}
enum DetectorMaskBit : std::int32_t {
  kDetPrimaryVertex = 1 << 1,   ///< bit 2, primary vertex used in the fit
  kDetVd            = 1 << 3,   ///< bit 4. **VD or VFT**, not VD alone
                                ///< (`mammoth.car:13556`); at LEP2 it is set on
                                ///< every VFT-reconstructed track, so a "has VD
                                ///< hits" cut written against it silently widens
  kDetId            = 1 << 4,   ///< bit 5
  kDetTpc           = 1 << 5,   ///< bit 6
  kDetRib           = 1 << 6,   ///< bit 7, barrel RICH
  kDetOd            = 1 << 8,   ///< bit 9
  kDetHpc           = 1 << 9,   ///< bit 10
  kDetHab           = 1 << 13,  ///< bit 14, barrel HCAL
  kDetMub           = 1 << 14,  ///< bit 15, barrel muon chambers
  kDetStic          = 1 << 18,  ///< bit 19
  kDetFca           = 1 << 20,  ///< bit 21
  kDetRif           = 1 << 21,  ///< bit 22, forward RICH
  kDetHaf           = 1 << 22,  ///< bit 23, forward HCAL
  kDetFcb           = 1 << 25,  ///< bit 26
  kDetEmf           = 1 << 26,  ///< bit 27
  kDetMuf           = 1 << 30   ///< bit 31, forward muon chambers
};

/// Which `edm4hep::TrackState::location` values appear on `TRAC_Tracks`, and
/// what each is here. Values are taken from EDM4hep, not restated, so this
/// cannot drift from the schema.
///
/// @collection{TRAC_Tracks}
enum TrackStateUse : std::int32_t {
  kStateAtOther       = edm4hep::TrackState::AtOther,
                        ///< PA.TRAX, every other surface: TOF(11), muon
                        ///< chambers (14, 17, 30)
  kStateAtIp          = edm4hep::TrackState::AtIP,
                        ///< the perigee helix; exactly one per track
  kStateAtFirstHit    = edm4hep::TrackState::AtFirstHit,
                        ///< PA.TRAX, TANAGRA detector id 0: the track's own
                        ///< first measured point
  kStateAtCalorimeter = edm4hep::TrackState::AtCalorimeter,
                        ///< PA.TRAX calorimeter crossings: HPC(9), HAB(13),
                        ///< HAF(22), EMF(26)
  kStateAtVertex      = edm4hep::TrackState::AtVertex
                        ///< AABTAG's impact parameter at AABTAG's own primary
                        ///< vertex; absent on tracks AABTAG skipped
};

/// Bits of the `MAIN_Particles` selection flag (SKELANA's `LVLOCK`), as masks.
///
/// Bit 32 is the sign bit, so a track locked by IFLVEC reads as a large
/// negative number: `flag > 0` and a bare `if (flag)` both mislead. Test
/// against these masks on the value cast to `std::uint32_t`. `-1` is not a bit
/// pattern at all but the converter's sentinel for a particle with no VECP
/// entry.
///
/// @collection{VECP_Particles_SelectionFlag}
enum SelectionFlagBit : std::uint32_t {
  kSelRejected = 0x00000001u,  ///< bit 1; failed the track-selection cuts,
                               ///< which `IFLSTR = 2` keeps and flags
  kSelIflvec   = 0x80000000u   ///< bit 32; locked by IFLVEC's class choice --
                               ///< ours is 22, "lock the charged outgoing"
};

class TrackingWriter : public CollectionWriter {
public:
  using CollectionWriter::CollectionWriter;
  void emit() override;
};

}  // namespace delphi_edm4hep::tracking
