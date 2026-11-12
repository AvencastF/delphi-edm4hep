// Truth domain — TruthGenWriter + TruthRecoLinkWriter.
//
// PSCLUJ (filled by PSHLUJ on sDST or PSFLUJ on fDST per bank-presence
// gating) feeds TruthGenWriter, with the DELSIM vertex record converted from
// cm to mm and a vertex/endpoint status per entry. PSCTBL exact tables (IPAST
// then ISTLU) feed TruthRecoLinkWriter, replacing the legacy helix-NN match.

#include "delphi_edm4hep/Truth/Truth.h"

#include "phdst/uxcom.hpp"
#include "phdst/uxlink.hpp"
#include "skelana/pscluj.hpp"
#include "skelana/psctbl.hpp"

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>
#include <podio/UserDataCollection.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

// JETSET helpers. Charges via LUCHGE in units of e/3 (exhaustive PDG
// coverage; replaces the 16-entry whitelist of the legacy code).
extern "C" {
  void pshluj_();
  void psfluj_();
  int  luchge_(int* kf);
}

namespace ph = phdst;
namespace sk = skelana;

namespace delphi_edm4hep::truth {

namespace {

inline float charge_from_pdg(int pdg) {
  return static_cast<float>(luchge_(&pdg)) / 3.0f;
}

// Gate the PSCLUJ unpacker on bank presence in the LDTOP chain, invoke
// the sDST or fDST routine, return sk::NP. (full-DST SH chain at
// LDTOP-3, sDST condensed LU at LDTOP-28/-29.)
int unpack_lujets() {
  sk::NP = 0;
  if (ph::LDTOP <= 0) return 0;
  const bool short_sim = (ph::IQ(ph::LDTOP - 2) > 28)
                      && (ph::LQ(ph::LDTOP - 28) != 0)
                      && (ph::LQ(ph::LDTOP - 29) != 0);
  const bool full_sim  = (ph::LQ(ph::LDTOP - 3) != 0);
  if      (short_sim) pshluj_();
  else if (full_sim)  psfluj_();
  return sk::NP;
}

// The PSCLUJ vertex VP(i,1..3) is a DELSIM vertex-bank (PV) position in cm,
// already in the DELSIM/reco frame (the primary sits at the simulated beam
// spot). PSFLUJ copies Q(LPV+5..7) unconverted (skelana.car PSFLUJ, 7122-7132);
// DELPHI's own dstana10 SPATAU divides differences of the same words by c in
// cm/ps. The PSCLUJ header's "mm" is the LUJETS convention, not this record's.
constexpr double kCm2Mm = 10.0;
constexpr double kNaN   = std::numeric_limits<double>::quiet_NaN();

// sDST_LUJ_GenParticles_VertexStatus / _EndpointStatus values.
enum VertexStatus : std::int32_t {
  kVtxOwn       = 0,   // its own DELSIM vertex (SH -> ST -> PV), VP x 10
  kVtxInherited = 1,   // no ST: DELSIM's point for it, its mother's decay point
  kVtxPrimary   = 2,   // no ST and no mother: the simulated primary vertex
  kVtxUnknown   = 3,   // not determined: NaN
};
enum EndpointStatus : std::int32_t {
  kEndDecay      = 0,  // the DELSIM decay vertex its own-vertex descendants start at
  kEndZeroFlight = 1,  // no ST: DELSIM decays it where it is produced
  kEndNone       = 2,  // no daughters in the record: NaN
  kEndUnknown    = 3,  // daughters, but no consistent decay vertex: NaN
};

// Per LU index (0-based), whether PSFLUJ/PSHLUJ filled VP from a vertex bank
// of the particle's own: fullDST, the SH has an ST whose PV exists (PSFLUJ
// 7122-7132, walked in the same SH order); shortDST, the SH is one of the
// SH-with-ST entries and its PV index is non-zero (PSHLUJ 5374-5398). For
// every other entry SKELANA leaves VP as PSBEG zeroed it (skelana.car 629) or
// sets it to 0 itself: "no vertex", not the origin. Gated as unpack_lujets().
std::vector<char> own_sim_vertex(int nGen) {
  std::vector<char> own(static_cast<std::size_t>(std::max(nGen, 0)), 0);
  if (nGen <= 0 || ph::LDTOP <= 0) return own;
  const bool short_sim = (ph::IQ(ph::LDTOP - 2) > 28)
                      && (ph::LQ(ph::LDTOP - 28) != 0)
                      && (ph::LQ(ph::LDTOP - 29) != 0);
  if (short_sim) {
    const int lstsh = ph::LQ(ph::LDTOP - 29);
    if (lstsh <= 0) return own;
    const int nsh0   = ph::IQ(lstsh + 2) / 1000;
    const int nstsh  = ph::IQ(lstsh + 3) / 1000;
    const int nst0   = ph::IQ(lstsh + 4) / 1000;
    const int lenhed = ph::IQ(lstsh + 1) - (nsh0 + 2 * nstsh + nst0) + 1;
    for (int i = 1; i <= nGen && i <= sk::NVMAX; ++i) {
      const int j = sk::ILUSH(i);
      if (j <= nsh0 || j > nsh0 + nstsh) continue;
      const int word = ph::IQ(lstsh + lenhed + nsh0 + 2 * (j - nsh0));
      own[static_cast<std::size_t>(i - 1)] = ((word >> 8) & 0xFF) != 0;  // JBYT(.,9,8)
    }
    return own;
  }
  int i = 0;
  for (int lsh = ph::LQ(ph::LDTOP - 3); lsh > 0 && i < nGen; lsh = ph::LQ(lsh)) {
    if (std::lround(ph::Q(lsh + 1)) == 0) continue;
    const int lst = ph::LQ(lsh - 4);
    own[static_cast<std::size_t>(i++)] = lst > 0 && ph::LQ(lst + 1) > 0;
  }
  return own;
}

// Per-event MC truth interaction point (sim-PV), in mm. Handles BOTH sim
// layouts: shortDST (LPVS at LDTOP-28, first sim-PV x/y/z at Q(ip+4..6) cm) and
// fullDST (LDTOP-3 -> LSH -> LST -> LPV, x/y/z at Q(lpv+5..7) cm); same decode
// as delphi-raw-nanoaod fillSimPV. Returns nullopt when no sim banks are present
// (real data / missing bank). It is the production vertex of truth entries that
// have neither a vertex of their own nor a mother (beams, Z/gamma*).
std::optional<std::array<double, 3>> read_sim_pv_mm() {
  if (ph::LDTOP <= 0) return std::nullopt;
  // shortDST sim layout: LPVS at LDTOP-28, x/y/z at +4/+5/+6 (cm). Guard the
  // down-link index exactly as unpack_lujets() does -- LDTOP-28 is only a valid
  // link when the LDTOP bank has >28 down-links -- so this is safe on real data
  // (and on fullDST, where it falls through to the LSH-chain walk below).
  if (ph::IQ(ph::LDTOP - 2) > 28) {
    const int lpvs = ph::LQ(ph::LDTOP - 28);
    if (lpvs > 0) {
      const int npvs = ph::IQ(lpvs + 1);
      if (npvs >= 1) {
        const int ip = lpvs + 1 + npvs;
        return std::array<double, 3>{ph::Q(ip + 4) * kCm2Mm,   // cm -> mm
                                     ph::Q(ip + 5) * kCm2Mm,
                                     ph::Q(ip + 6) * kCm2Mm};
      }
    }
  }
  // fullDST sim layout (our DELSIM SDST takes this path): walk the LSH chain at
  // LDTOP-3, first SH with Q(lsh+1)!=0 -> LST at LSH-4 -> LPV at LST+1; x/y/z at
  // LPV+5/+6/+7 (cm). Same decode as raw-nanoaod fillSimPV (PSFLUJ, skelana.car ~L7088).
  for (int lsh = ph::LQ(ph::LDTOP - 3); lsh > 0; lsh = ph::LQ(lsh)) {
    if (std::lround(ph::Q(lsh + 1)) == 0) continue;
    const int lst = ph::LQ(lsh - 4);
    if (lst <= 0) continue;
    const int lpv = ph::LQ(lst + 1);
    if (lpv <= 0) continue;
    return std::array<double, 3>{ph::Q(lpv + 5) * kCm2Mm,   // cm -> mm
                                 ph::Q(lpv + 6) * kCm2Mm,
                                 ph::Q(lpv + 7) * kCm2Mm};
  }
  return std::nullopt;
}

}  // namespace

// ---------------------------------------------------------------------------
void TruthGenWriter::emit() {
  GenParticleResult result;

  const int nGen = unpack_lujets();
  const std::size_t n = static_cast<std::size_t>(std::max(nGen, 0));
  edm4hep::MCParticleCollection mc;
  podio::UserDataCollection<std::int32_t> vertexStatusCol, endpointStatusCol;
  result.handles.reserve(n);

  // Vertices and endpoints follow DELSIM's own placement (delsim36.car
  // SXSHST/SXSDK): a particle with an ST bank starts at that ST's vertex; a
  // decaying one (SXSDK) ends at a secondary vertex from which the STs of its
  // final-state descendants hang -- including those reached through daughters
  // without an ST (e.g. a generator-decayed pi0 that DELSIM does not decay
  // itself has no ST; its photons start at the tau's decay vertex). A
  // particle without an ST is therefore produced at its mother's decay point
  // and decays where it is produced; with no mother either, it is at the
  // simulated primary. Real data (nGen = 0) never reads sim banks -> empty
  // truth.
  using V3 = std::array<double, 3>;
  const V3 nan3{kNaN, kNaN, kNaN};
  const auto own = own_sim_vertex(nGen);
  const auto spv = nGen >= 1 ? read_sim_pv_mm() : std::nullopt;
  auto mother = [&](int i) {   // 0-based mother index, or -1
    const int m = sk::KP(i + 1, 3);
    return (m >= 1 && m <= nGen && m != i + 1) ? m - 1 : -1;
  };
  // Nearest own-vertex ancestor of i (0-based), or -1; at most n steps, so a
  // malformed mother cycle cannot hang the walk.
  auto own_ancestor = [&](int i) {
    int a = mother(i);
    for (int steps = 0; a >= 0 && !own[static_cast<std::size_t>(a)]; ++steps) {
      if (steps > nGen) return -2;
      a = mother(a);
    }
    return a;
  };
  std::vector<V3> vtx(n, nan3), end(n, nan3);
  std::vector<std::int32_t> vstat(n, kVtxUnknown), estat(n, kEndNone);
  std::vector<char> has_daughter(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const int m = mother(static_cast<int>(i));
    if (m >= 0) has_daughter[static_cast<std::size_t>(m)] = 1;
    if (own[i]) {
      const int lu = static_cast<int>(i) + 1;
      vtx[i] = {sk::VP(lu, 1) * kCm2Mm, sk::VP(lu, 2) * kCm2Mm, sk::VP(lu, 3) * kCm2Mm};
      vstat[i] = kVtxOwn;
    }
  }
  // Decay vertex of each own-vertex particle: the start of every own-vertex
  // descendant whose nearest own-vertex ancestor it is. They hang from one
  // DELSIM vertex; any disagreement leaves the endpoint undetermined.
  for (std::size_t d = 0; d < n; ++d) {
    if (!own[d]) continue;
    const int a = own_ancestor(static_cast<int>(d));
    if (a < 0) continue;
    const auto ai = static_cast<std::size_t>(a);
    if (estat[ai] == kEndNone) { end[ai] = vtx[d]; estat[ai] = kEndDecay; }
    else if (estat[ai] == kEndDecay && end[ai] != vtx[d]) { end[ai] = nan3; estat[ai] = kEndUnknown; }
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (own[i]) {
      if (has_daughter[i] && estat[i] == kEndNone) estat[i] = kEndUnknown;  // NaN
      continue;
    }
    // Production point of an entry without an ST: the decay vertex of its
    // nearest own-vertex ancestor, else (no such ancestor) the primary.
    const int a = own_ancestor(static_cast<int>(i));
    if (a >= 0) {
      vtx[i] = end[static_cast<std::size_t>(a)];
      vstat[i] = std::isnan(vtx[i][0]) ? kVtxUnknown : kVtxInherited;
    } else if (a == -1 && spv) {
      vtx[i] = *spv;
      vstat[i] = mother(static_cast<int>(i)) < 0 ? kVtxPrimary : kVtxInherited;
    }
    if (has_daughter[i]) {
      end[i] = vtx[i];
      estat[i] = vstat[i] == kVtxUnknown ? kEndUnknown : kEndZeroFlight;
    }
  }

  // KP(i,1)=status, KP(i,2)=PDG, PP(i,1..3)=p, PP(i,5)=mass (GeV). No time
  // is read by SKELANA (VP(i,4) = 0, unfilled): published as NaN.
  for (std::size_t i = 0; i < n; ++i) {
    const int lu = static_cast<int>(i) + 1;
    auto mp = mc.create();
    mp.setPDG(sk::KP(lu, 2));
    mp.setGeneratorStatus(static_cast<std::int16_t>(sk::KP(lu, 1)));
    mp.setMomentum({sk::PP(lu, 1), sk::PP(lu, 2), sk::PP(lu, 3)});
    mp.setMass(sk::PP(lu, 5));
    mp.setVertex({vtx[i][0], vtx[i][1], vtx[i][2]});
    mp.setEndpoint({end[i][0], end[i][1], end[i][2]});
    mp.setTime(static_cast<float>(kNaN));
    mp.setCharge(charge_from_pdg(sk::KP(lu, 2)));
    vertexStatusCol.push_back(vstat[i]);
    endpointStatusCol.push_back(estat[i]);
    result.handles.push_back(mp);
  }

  // Parent/daughter graph.
  for (std::size_t i = 0; i < n; ++i) {
    const int m = mother(static_cast<int>(i));
    if (m < 0) continue;
    result.handles[i].addToParents(result.handles[static_cast<std::size_t>(m)]);
    result.handles[static_cast<std::size_t>(m)].addToDaughters(result.handles[i]);
  }

  put(std::move(mc), "LUJ", "GenParticles", Provenance::Derived);
  put(std::move(vertexStatusCol),   "LUJ", "GenParticles_VertexStatus",   Provenance::Derived);
  put(std::move(endpointStatusCol), "LUJ", "GenParticles_EndpointStatus", Provenance::Derived);
  ctx_.gen_truth = std::move(result);
}

// ---------------------------------------------------------------------------
void TruthRecoLinkWriter::emit() {
  edm4hep::RecoMCParticleLinkCollection links;

  // Need both upstream writers' outputs.
  if (!ctx_.gen_truth || !ctx_.tracking) {
    put(std::move(links), "TBL", "RecoToGen", Provenance::Transcribed);   // emit empty + return
    return;
  }
  const auto& gen      = *ctx_.gen_truth;
  const auto& tracking = *ctx_.tracking;

  // PSCTBL.NPA = # PA particles = VECP entries. For each VECP index j:
  //   ist = IPAST(j)                 (ST index; 0 if no MC ancestor)
  //   ilu = ISTLU(ist)               (LU index)
  //   particle_idx = vecp_to_particle[j]   (-1 if Tracking dropped it)
  const int nPA   = sk::NPA();
  const int nGen  = static_cast<int>(gen.handles.size());
  const auto& v2p = tracking.vecp_to_particle;

  for (int j = 1; j <= nPA; ++j) {
    if (j >= static_cast<int>(v2p.size())) break;
    const int particle_idx = v2p[j];
    if (particle_idx < 0) continue;
    const int ist = sk::IPAST(j);
    if (ist <= 0) continue;
    const int ilu = sk::ISTLU(ist);
    if (ilu <= 0 || ilu > nGen) continue;

    auto link = links.create();
    link.setFrom(tracking.particle_handles[particle_idx]);
    link.setTo  (gen.handles[ilu - 1]);
    link.setWeight(1.0f);
  }

  put(std::move(links), "TBL", "RecoToGen", Provenance::Transcribed);
}

}  // namespace delphi_edm4hep::truth
