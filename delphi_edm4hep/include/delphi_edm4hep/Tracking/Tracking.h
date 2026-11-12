// Tracking domain: PA.MAIN + PA.TRAC + PSCVEC (VECP / LVLOCK).
// TrackingWriter emits:
//   <prefix>_TRAC_Tracks                   Track + AtIP TrackState + 5x5 cov
//   <prefix>_MAIN_Particles                charged + neutral, VECP 4-mom;
//                                    NaN 4-mom for a charged PA without a
//                                    charged VECP row, plus SKELANA's
//                                    recovered neutral as its own entry
//   <prefix>_VECP_Particles_SelectionFlag  UserData<int32>
//   <prefix>_MAIN_Particles_VecpStatus     UserData<int32>, VecpStatus
//                                    (TrackingData.h)
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
