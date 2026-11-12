// Calorimeter domain: EMNC (label 22) + HCNC (label 23) / HCAL (label 3) on sDST.
// CalorimeterWriter reads ctx_.tracking and emits:
//   <tag>_EMNC_Showers   ECAL (HPC bit 0 / EMF bit 1 in Cluster.type)
//   <tag>_HCNC_Showers   HCAL (bit 2), HACCOR-corrected
//   <tag>_HCAL_Showers   HCAL (bit 2), as reconstructed
// Sets Particle.clusters on each owner Particle via
// tracking.calorimeterOwner(paIdx): the PA's own Particle, or the
// SKELANA-recovered neutral built from a charged PA's showers.

#pragma once

#include "delphi_edm4hep/CollectionWriter.h"

namespace delphi_edm4hep::calorimeter {

class CalorimeterWriter : public CollectionWriter {
public:
  using CollectionWriter::CollectionWriter;
  void emit() override;
};

}  // namespace delphi_edm4hep::calorimeter
