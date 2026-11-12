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
//                                    one
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

namespace delphi_edm4hep::tracking {

class TrackingWriter : public CollectionWriter {
public:
  using CollectionWriter::CollectionWriter;
  void emit() override;
};

}  // namespace delphi_edm4hep::tracking
