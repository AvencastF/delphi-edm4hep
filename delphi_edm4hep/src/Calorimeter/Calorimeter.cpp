// Calorimeter domain — shower clusters from the PA calorimeter modules.
//
// Walks the PA chain and emits one Cluster per shower. Cluster.type bits
// encode the sub-detector:
//   bit 0 = HPC   (EMNC, idet 9)
//   bit 1 = FEMC  (EMNC, idet 26)
//   bit 2 = HCAL  (HCNC or HCAL module)
//
// Hadron showers come from two modules, and which the production wrote is a
// per-processing choice: HCNC(23), HCAL(3), or both. Both are emitted. Where
// both are present HCNC supersedes HCAL — it is HCAL re-associated by HACCOR
// (ecorrxx.car:273, "HCNC module overwrites module 3 info"). The README
// states the consequence for analysis.
//
// One collection per source module. Charged/neutral is not split here — it
// is inferable from the owning Particle's `clusters` relation.
//
// Module layouts are documented in the ShortDST content note, DELPHI 97-146
// PROG-221, doi:10.7483/OPENDATA.DELPHI.LJAQ.LMZ0, under "PA extra-module
// EMNC(22)" and "HCNC(23)". Word numbers below refer to those tables.

#include "delphi_edm4hep/Calorimeter/Calorimeter.h"

#include "delphi_edm4hep/internal/BankReader.h"

#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/MutableCluster.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <podio/Frame.h>

#include <algorithm>
#include <cstdint>
#include <string>

namespace delphi_edm4hep::calorimeter {

namespace {

// DELPHI cm -> EDM4hep mm.
constexpr double kCm2Mm = 10.0;

// Cluster.type bit encoding (see file header).
constexpr std::int32_t kTypeBitHPC  = 0x01;
constexpr std::int32_t kTypeBitEMF  = 0x02;
constexpr std::int32_t kTypeBitHCAL = 0x04;

// Detector identifiers carried in EMNC word +2.
constexpr int kIdetHPC = 9;
constexpr int kIdetEMF = 26;

// Defensive bound on the per-shower layer count. Both modules document
// layers 1..4.
constexpr int kMaxHcalLayers = 200;

// Shower header: energy (GeV) at word `e`, then x, y, z (cm).
void setHeaderEnergyPosition(edm4hep::MutableCluster clu,
                             const banks::Record& shower, int e) {
  clu.setEnergy(shower.real(e));
  clu.setPosition({
    static_cast<float>(shower.real(e + 1) * kCm2Mm),
    static_cast<float>(shower.real(e + 2) * kCm2Mm),
    static_cast<float>(shower.real(e + 3) * kCm2Mm),
  });
}

// Electromagnetic showers for one PA. Emits one Cluster per shower and
// attaches each to the owning particle.
void walkEMNC(int lpa,
              edm4hep::MutableReconstructedParticle pfo,
              edm4hep::ClusterCollection& cluCol,
              const std::vector<std::vector<edm4hep::CalorimeterHit>>* pa_hits)
{
  const auto emnc = banks::find(lpa, "EMNC");
  if (!emnc) return;

  // Module word +2 packs both counts: nshowr + 100 * idet. idet is uniform
  // within one module.
  const int nsidet = emnc->integer(2);
  const int nshowr = nsidet % 100;
  const int idet   = nsidet / 100;
  if (nshowr <= 0) return;

  const std::int32_t type_bit =
      (idet == kIdetHPC) ? kTypeBitHPC :
      (idet == kIdetEMF) ? kTypeBitEMF :
                           0;          // unknown detector — emit anyway

  // The HPC has ten layers.
  constexpr int kMaxHpcLayers = 10;

  auto shower = emnc->subRecord(2);
  for (int ns = 0; ns < nshowr; ++ns) {
    auto clu = cluCol.create();
    clu.setType(type_bit);

    // Calorimeter hits of the same shower, decoded by EmcaWriter from the
    // same PA. PA.EMCA and PA.EMNC list a PA's showers in the same order.
    if (pa_hits && ns < static_cast<int>(pa_hits->size())) {
      for (const auto& hit : (*pa_hits)[static_cast<std::size_t>(ns)]) {
        clu.addToHits(hit);
      }
    }
    setHeaderEnergyPosition(clu, shower, 1);

    int stride = 4;               // FEMC carries no layer detail here

    if (idet == kIdetHPC) {
      // HPC continuation: +5 = layers hit, +6 = layer bit pattern,
      // +7.. = energy of each layer present, in pattern order.
      const int nlay = shower.integer(5);
      const int patt = shower.integer(6);
      if (nlay > 0 && nlay <= kMaxHpcLayers) {
        // Expanded to a dense ten-element profile, zero for layers absent
        // from the pattern, so subdetectorEnergies has a fixed length per
        // HPC cluster.
        int filled = 0;
        for (int nl = 1; nl <= kMaxHpcLayers && filled < nlay; ++nl) {
          const bool hit = ((patt >> (nl - 1)) & 1) != 0;
          float layerE = 0.f;
          if (hit) {
            ++filled;
            layerE = shower.real(6 + filled);
          }
          clu.addToSubdetectorEnergies(layerE);
        }
      }
      stride = 6 + nlay;
    }
    // FEMC per-block detail lives in PA module EL(5), which is fullDST only.

    pfo.addToClusters(clu);
    shower = shower.subRecord(stride);
  }
}

// Hadronic showers for one PA. Per-layer energies go to
// subdetectorEnergies and the matching layer indices to shapeParameters,
// so the two stay index-parallel without a second collection.
void walkHCNC(int lpa,
              edm4hep::MutableReconstructedParticle pfo,
              edm4hep::ClusterCollection& cluCol)
{
  const auto hcnc = banks::find(lpa, "HCNC");
  if (!hcnc) return;

  // Module word +2 = number of showers.
  const int nshowr = hcnc->integer(2);
  if (nshowr <= 0) return;

  auto shower = hcnc->subRecord(2);
  for (int ns = 0; ns < nshowr; ++ns) {
    auto clu = cluCol.create();
    clu.setType(kTypeBitHCAL);
    setHeaderEnergyPosition(clu, shower, 1);

    // +5 = layer hits. Each hit is two words:
    //   +5 + 2n - 1   energy
    //   +5 + 2n       1000 * layer + channels in that layer
    const int nlay = shower.integer(5);
    if (nlay > 0 && nlay < kMaxHcalLayers) {
      for (int nl = 1; nl <= nlay; ++nl) {
        const float hitE  = shower.real(5 + 2 * nl - 1);
        const int   layer = shower.integer(5 + 2 * nl) / 1000;
        clu.addToSubdetectorEnergies(hitE);
        clu.addToShapeParameters(static_cast<float>(layer));
      }
    }

    pfo.addToClusters(clu);
    shower = shower.subRecord(5 + 2 * std::max(0, nlay));
  }
}

// Hadronic showers from the HCAL module. Same content as walkHCNC at
// different offsets, so the emitted fields match. Showers are named
// sub-records rather than a fixed stride (PSHHAC, skelana.car:3504).
void walkHCAL(int lpa,
              edm4hep::MutableReconstructedParticle pfo,
              edm4hep::ClusterCollection& cluCol)
{
  const auto hcal = banks::find(lpa, "HCAL");
  if (!hcal) return;

  // Module word +2 = number of showers.
  const int nshowr = hcal->integer(2);

  for (int ns = 1; ns <= nshowr; ++ns) {
    const auto shower = hcal->subRecord("HCAL.SHOWER", ns);
    if (!shower) break;

    auto clu = cluCol.create();
    clu.setType(kTypeBitHCAL);
    setHeaderEnergyPosition(clu, *shower, 3);

    // +10 = layers hit, +11 = layer bit pattern. Each hit is two words:
    //   +11 + 2n - 1   energy
    //   +11 + 2n       1000 * layer + active channels in that layer
    const int nlay = shower->integer(10);
    if (nlay > 0 && nlay < kMaxHcalLayers) {
      for (int nl = 1; nl <= nlay; ++nl) {
        clu.addToSubdetectorEnergies(shower->real(11 + 2 * nl - 1));
        clu.addToShapeParameters(
            static_cast<float>(shower->integer(11 + 2 * nl) / 1000));
      }
    }

    pfo.addToClusters(clu);
  }
}

}  // namespace

void CalorimeterWriter::emit()
{
  edm4hep::ClusterCollection emncCol;
  edm4hep::ClusterCollection hcncCol;
  edm4hep::ClusterCollection hcalCol;

  if (ctx_.tracking) {
    const auto& tracking = *ctx_.tracking;
    pawalk::forEachPA([&](int lpa, int paIdx) {
      // Usually the PA's own particle; for a charged PA that SKELANA
      // recovered as a neutral, the recovered neutral made from these
      // showers, so each shower stays attached to exactly one particle.
      const int particle_idx = tracking.calorimeterOwner(paIdx);
      if (particle_idx < 0) return;
      auto pfo = tracking.particle_handles[particle_idx];

      const std::vector<std::vector<edm4hep::CalorimeterHit>>* pa_hits = nullptr;
      if (ctx_.emca && paIdx < static_cast<int>(ctx_.emca->pa_shower_hits.size())) {
        pa_hits = &ctx_.emca->pa_shower_hits[static_cast<std::size_t>(paIdx)];
      }
      walkEMNC(lpa, pfo, emncCol, pa_hits);
      walkHCNC(lpa, pfo, hcncCol);
      walkHCAL(lpa, pfo, hcalCol);
    });
  }

  // Emitted even when empty, so the schema does not depend on whether
  // tracking ran.
  put(std::move(emncCol), "EMNC", "Showers", Provenance::Transcribed);
  put(std::move(hcncCol), "HCNC", "Showers", Provenance::Transcribed);
  put(std::move(hcalCol), "HCAL", "Showers", Provenance::Transcribed);
}

}  // namespace delphi_edm4hep::calorimeter
