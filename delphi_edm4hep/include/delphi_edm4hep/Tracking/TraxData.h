// TraxData.h — value type passed through EventContext.
// Separate from Trax.h so CollectionWriter.h can include this without pulling
// in the writer class.

#pragma once

#include <edm4hep/TrackState.h>

#include <vector>

namespace delphi_edm4hep::trax {

// Output of TraxWriter::emit(), indexed by PA-walk index: the extrapolation
// states of each PA's track. The track writers append them to the track they
// build for that PA, while it is still mutable.
struct Output {
  std::vector<std::vector<edm4hep::TrackState>> pa_to_states;
};

}  // namespace delphi_edm4hep::trax
