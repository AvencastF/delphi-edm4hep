// Tracking domain — implementation S3 (pass-1).
//
// Walks the PA chain (LDTOP-1 -> per-PV -> per-PA), emitting:
//   <prefix>_TRAC_Tracks    (Track + TrackState[AtIP] + 5x5 helix-basis cov)
//   <prefix>_MAIN_Particles (charged + neutral; 4-mom from sk::VECP)
//   <prefix>_VECP_Particles_SelectionFlag       (UserData int32)
//   <prefix>_VECP_Particles_SelectionWord       (UserData int32, LVSELE)
//   <prefix>_MAIN_Particles_VecpStatus         (UserData int32, VecpStatus)
//   <prefix>_MAIN_Particles_ChargeCode          (UserData int32, PA.MAIN charge code)
//   <prefix>_VECP_Particles_Charge              (UserData float, VECP(7))
//   <prefix>_MAIN_Particles_ReconstructionCode  (UserData int32)
//   <prefix>_MAIN_Particles_IdentificationWord  (UserData int32, IQ(LPA+3))
//   <prefix>_MAIN_Particles_DetectorMask        (UserData int32)
//   <prefix>_MAIN_Particles_MassCode            (UserData int32, PA.MAIN word 10)
//   <prefix>_MAIN_Particles_EnergyError         (UserData float, PA.MAIN word 19)
//
// Each Track links to the track elements from its PA (see TrackElements).
//   <prefix>_MAIN_Particles_TrackLength         (UserData float, cm)
//   <prefix>_QTRAC_Tracks_d0PV / _z0PV / _d0BS  (UserData float; QTRAC 38..40)
//   <prefix>_QTRAC_Tracks_Perigee*              (UserData float; QTRAC 4..23)
//   <prefix>_QTRAC_Tracks_PerigeeRefitted       (UserData int32)
//   <prefix>_QTRAC_Tracks_FitQuality            (UserData float; QTRAC 28)
//   <prefix>_QTRAC_Tracks_EnergyError           (UserData float; QTRAC 41)
//
// No shape moments, no sigma calibration, no perigee-momentum fallback —
// bank-truth values only (those custom operations are intentionally
// dropped). A value with no source is NaN, never 0.

#include "delphi_edm4hep/Tracking/Tracking.h"

#include "delphi_edm4hep/internal/AabtagTrackState.h"

#include "skelana/pscbsp.hpp"
#include "skelana/pscflg.hpp"

#include "delphi_edm4hep/Helix.h"
#include "delphi_edm4hep/internal/PaWalk.h"
#include "skelana/pscvec.hpp"
#include "skelana/psctra.hpp"
#include "skelana/pscvtx.hpp"
#include "skelana/trackntp.hpp"

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <podio/Frame.h>
#include <podio/UserDataCollection.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <string>

namespace ph = phdst;
namespace sk = skelana;
namespace aa = delphi_edm4hep::aabtag;

namespace {
constexpr double kCm2Mm = 10.0;
constexpr float  kNotMeasured = std::numeric_limits<float>::quiet_NaN();

// QTRAC(4..23), SKELANA's track parameters in PA.TRAC words 2..21 layout:
// perigee (d0, z0, theta, phi, 1/R) then the 15-word lower-triangular weight
// matrix of those five, row by row. Published raw, in DELPHI units.
constexpr int kPerigeeWords = 20;
constexpr std::array<const char*, kPerigeeWords> kPerigeeNames = {
    "Tracks_PerigeeD0",       "Tracks_PerigeeZ0",       "Tracks_PerigeeTheta",
    "Tracks_PerigeePhi",      "Tracks_PerigeeInvR",
    "Tracks_PerigeeWeight11", "Tracks_PerigeeWeight21", "Tracks_PerigeeWeight22",
    "Tracks_PerigeeWeight31", "Tracks_PerigeeWeight32", "Tracks_PerigeeWeight33",
    "Tracks_PerigeeWeight41", "Tracks_PerigeeWeight42", "Tracks_PerigeeWeight43",
    "Tracks_PerigeeWeight44", "Tracks_PerigeeWeight51", "Tracks_PerigeeWeight52",
    "Tracks_PerigeeWeight53", "Tracks_PerigeeWeight54", "Tracks_PerigeeWeight55"};
}  // namespace

namespace delphi_edm4hep::tracking {

namespace {

// The DELPHI perigee -> EDM4hep helix conversion (Jacobian push-forward,
// weight-matrix inversion, TrackState emit) lives in the shared
// delphi_edm4hep::Helix class (Helix.h) so PA.TRAC and PA.ELTR use one
// definition.

// Look up the SKELANA VECP index for a given PA's lpa, using the LVECP
// back-link. The `want_charged` flag gates: a charged-track PA's own row
// has VECP(7,*) != 0; a neutral row has VECP(7,*) == 0. A neutral row
// whose LVECP is a charged PA is the one PSHCHREC re-enters after
// calorimeter recovery (the only other source of extra neutral rows,
// REMCLU, is off: IFLECL = 22). Returns 0 if no match.
int find_vecp_index(int lpa, bool want_charged)
{
  for (int i = sk::LVPART; i <= sk::NVECP; ++i) {
    if (sk::LVECP(i) != lpa) continue;
    const int q = static_cast<int>(std::lround(sk::VECP(7, i)));
    if (want_charged && q != 0) return i;
    if (!want_charged && q == 0) return i;
  }
  return 0;
}

}  // namespace

void TrackingWriter::emit()
{
  // Storage that will be moved into the Frame at the end. Handles taken
  // before the move remain valid (podio guarantees handle stability
  // across the collection move).
  edm4hep::TrackCollection                 trkCol;
  edm4hep::ReconstructedParticleCollection pfoCol;
  podio::UserDataCollection<std::int32_t>  lvlockCol;
  podio::UserDataCollection<std::int32_t>  lvseleCol;
  podio::UserDataCollection<std::int32_t>  statusCol;
  podio::UserDataCollection<std::int32_t>  chargeCodeCol;
  podio::UserDataCollection<float>         vecpChargeCol;
  podio::UserDataCollection<std::int32_t>  codeCol;
  podio::UserDataCollection<std::int32_t>  idWordCol;
  podio::UserDataCollection<std::int32_t>  detCol;
  podio::UserDataCollection<std::int32_t>  massCodeCol;
  podio::UserDataCollection<float>         energyErrorCol;
  podio::UserDataCollection<float>         lengthCol;
  podio::UserDataCollection<float>         d0PvCol;
  podio::UserDataCollection<float>         z0PvCol;
  podio::UserDataCollection<float>         d0BsCol;
  std::array<podio::UserDataCollection<float>, kPerigeeWords> perigeeCols;
  podio::UserDataCollection<std::int32_t>  refittedCol;
  podio::UserDataCollection<float>         fitQualityCol;
  podio::UserDataCollection<float>         qtracEnergyErrorCol;

  Output result;
  result.particle_handles.reserve(static_cast<std::size_t>(sk::NVECP + 1));
  // VECP-index map: size sk::NVECP + 1, all -1 initially. Index 0 unused.
  result.vecp_to_particle.assign(static_cast<std::size_t>(sk::NVECP + 1), -1);
  // PA-index map: sized lazily as we walk; -1 entries for PAs we skip.
  // Reserve a reasonable upper bound to avoid mid-walk reallocation.
  result.pa_to_particle.reserve(static_cast<std::size_t>(sk::NVECP + 16));

  // Helper: register a freshly-created particle handle in the output
  // metadata. vecp_i = 0 means "no VECP match" (the handle still goes
  // into particle_handles, but vecp_to_particle is not updated).
  // The current paIdx is taken from the surrounding forEachPA scope.
  auto record_particle = [&](edm4hep::MutableReconstructedParticle pfo,
                              int vecp_i, int paIdx) {
    const int new_idx = static_cast<int>(result.particle_handles.size());
    result.particle_handles.push_back(pfo);
    if (vecp_i >= 1 && vecp_i < static_cast<int>(result.vecp_to_particle.size())) {
      result.vecp_to_particle[vecp_i] = new_idx;
    }
    if (paIdx >= static_cast<int>(result.pa_to_particle.size())) {
      result.pa_to_particle.resize(static_cast<std::size_t>(paIdx + 1), -1);
    }
    result.pa_to_particle[paIdx] = new_idx;
  };

  // Helper: push the BS/PV-corrected impact-parameter triplet for a
  // VECP-matched track, or NaN sentinels if vecp_i is 0.
  //
  // CONVENTION / FOOTGUN (QTRAC_Tracks_d0PV / _z0PV / _d0BS): sk::QTRAC values
  // converted cm -> **mm**, DELPHI sign (rho = Vx sinphi - Vy cosphi), parallel
  // to TRAC_Tracks (charged-only). PSHTRA fills them for data and MC alike.
  // It pre-fills -999 cm and keeps it when there is no usable primary vertex
  // (KVTX(16,1) > 0, for d0PV/z0PV) or TBDCAE gives no answer; that sentinel
  // is published as NaN. For the geometric impact parameter use
  // PV_Tracks_d0PV (Vertex.cpp): also **mm** but LCIO sign, so d0 differs by
  // a factor -1 while z0 agrees -- do not mix them.
  static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  auto qtrac_mm = [](int k, int vecp_i) {
    const float cm = sk::QTRAC(k, vecp_i);
    return cm == -999.f ? kNaN : cm * 10.f;  // cm -> mm
  };
  //
  // Perigee words: QTRAC(4..23) as stored. PSHTRA copies PA.TRAC words 2..21
  // there; PSHCTRECPBS, the only other writer, replaces them with its
  // primary-vertex-constrained refit (the VECP 4-vector, QTRAC 38..40 and the
  // VECP charge are refilled from the same refit). Refitted = 1 when the
  // stored words differ bit for bit from the track's PA.TRAC words, 0 when
  // they are identical, -1 without an own VECP row (words NaN).
  //
  // Fit words QTRAC(28) and QTRAC(41), as stored, NaN without an own row.
  // PSHTRA copies PA.MAIN word 16 (chi2 of the fit without VD) and word 19
  // (energy error, GeV) there; a kept refit replaces them with CONFPV's
  // PCHI, which is the chi2 *probability* PROB(chi2, ndf) >= 0.01 and not a
  // chi2, and with the refit's momentum error PERR (GeV). Refitted says which.
  auto push_qtrac_or_nan = [&](int vecp_i, int ltrac) {
    fitQualityCol.push_back(vecp_i >= 1 ? sk::QTRAC(28, vecp_i) : kNaN);
    qtracEnergyErrorCol.push_back(vecp_i >= 1 ? sk::QTRAC(41, vecp_i) : kNaN);
    if (vecp_i >= 1) {
      d0PvCol.push_back(qtrac_mm(38, vecp_i));
      z0PvCol.push_back(qtrac_mm(39, vecp_i));
      d0BsCol.push_back(sk::IERRBS == 0 ? qtrac_mm(40, vecp_i) : kNaN);
    } else {
      d0PvCol.push_back(kNaN);
      z0PvCol.push_back(kNaN);
      d0BsCol.push_back(kNaN);
    }
    bool refitted = false;
    for (int w = 0; w < kPerigeeWords; ++w) {
      const float stored = vecp_i >= 1 ? sk::QTRAC(4 + w, vecp_i) : kNaN;
      perigeeCols[w].push_back(stored);
      refitted = refitted || (vecp_i >= 1 &&
          std::bit_cast<std::uint32_t>(stored) !=
          std::bit_cast<std::uint32_t>(ph::Q(ltrac + 2 + w)));
    }
    refittedCol.push_back(vecp_i >= 1 ? refitted : -1);
  };

  using namespace delphi_edm4hep::pawalk;

  // Helper: extend pa_to_particle with -1 so the map stays parallel to
  // the PA-walk index even when we return early for a skipped PA.
  auto mark_pa_skipped = [&](int paIdx) {
    if (paIdx >= static_cast<int>(result.pa_to_particle.size())) {
      result.pa_to_particle.resize(static_cast<std::size_t>(paIdx + 1), -1);
    }
  };

  // Charged PAs emitted without a 4-vector of their own (VecpStatus
  // RecoveredSource or NoRow), and charged PAs that unexpectedly have both a
  // charged and a neutral VECP row. Logged once per event below so the
  // degradation is loud.
  int n_recovered = 0;
  int n_no_row = 0;
  int n_dual_rows = 0;

  // Per-particle words running parallel to <tag>_MAIN_Particles.
  //
  //   LVLOCK  selection verdict; 0 = selected. Neutrals are selected too
  //           (PSHSNT, with the per-table calorimeter thresholds), so they
  //           carry a real verdict rather than a sentinel. -1 = no usable
  //           verdict: no VECP row, or the verdict relied on an unavailable
  //           beam spot. VecpStatus says which.
  //   LVSELE  SKELANA's per-row selection word (TRKINF), raw: bit 1 and 32
  //           are the ones copied into LVLOCK; 2-11 the failed cuts (ICTFAIL
  //           or INTFAIL shifted up by one), 20-22 the recovery step. -1
  //           without an own VECP row (no real word has every bit set).
  //           Recovered neutrals: see recovery_slot below.
  //   status  VecpStatus of the entry.
  //   charge code  PA.MAIN Q(LMAIN+8) as stored: 0 neutral, 1 positive,
  //           2 negative, 3 undefined. The particle charge is derived from it.
  //   VECP charge  VECP(7) of the entry's own row as SKELANA stored it, NaN
  //           without one (VecpStatus 1 and 3). PSHVEC writes +1/-1 for code
  //           1/2 and 999 for any other charged code; a beam-spot refit
  //           (PSHCTRECPBS) rewrites it from the charge-code word of its
  //           PA.MAIN copy passed through CONFPV (2 -> -1).
  //   code    PXPHOT reconstruction code, IQ(LPA+3) bits 19-25. SKELANA reads
  //           it to reject VD-only tracks (75 with z, 77 without z, VFT-only
  //           from ISVER 107) and ID+VD-only tracks without z (72).
  //   id word IQ(LPA+3) whole; `code` is its bits 19-25.
  //   mass code  PA.MAIN Q(LMAIN+10), SKELANA's VECP(8), an integer code.
  //   energy error  PA.MAIN Q(LMAIN+19), SKELANA's QTRAC(41) before any
  //           refit, as stored (SKELANA defines no unset value for it).
  //   length  track length in cm, the one selection quantity that cannot be
  //           reconstructed from the other emitted collections.
  // charge code, code, id word, mask, mass code, energy error and length are
  // PA words; a RecoveredNeutral entry repeats those of the PA it was built
  // from.
  //
  // Converter-side correction of a SKELANA ordering bug (PSHCHREC). Each
  // charged track recovered as a neutral gets its own work word: the k-th
  // recovery (PSHCTRECCAL, charged pass, PV->PA order) fills
  // LVSELE(LVWORK+k) and LVLOCK(LVWORK+k). With IFLODR = 1, PSHORTTRK
  // re-enters the recovered tracks only after the neutral pass, slot k as
  // VECP row NVECP-INREC+k, and PSHCHREC then copies
  // LVSELE/LVLOCK(LVWORK+INREC) -- the LAST recovery's -- into every one of
  // them. The converter publishes slot k for row NVECP-INREC+k instead, and
  // rebuilds that row's LVLOCK as SKELANA's final PSHORTTRK loop would have:
  // the work LVLOCK plus the word's bits 1 and 32. With one recovered
  // neutral, k = INREC and nothing changes. The mapping holds only if the
  // recovered rows are exactly the last INREC rows; if they are not, the
  // event's recovered neutrals get -1 (no usable word) and a warning.
  const int first_recovered_row = sk::NVECP - sk::INREC + 1;
  bool recovery_slots_ok = sk::IFLODR == 1 && sk::INREC >= 1 &&
                           sk::INREC <= sk::MTRACK;
  if (recovery_slots_ok) {
    for (int i = sk::LVPART; i <= sk::NVECP; ++i) {
      const int lm = std::lround(sk::VECP(7, i)) == 0 && sk::LVECP(i) > 0
                         ? lphpa("MAIN", sk::LVECP(i)) : 0;
      const bool recovered = lm > 0 && std::lround(ph::Q(lm + 8)) != 0;
      if (recovered != (i >= first_recovered_row)) recovery_slots_ok = false;
    }
  }
  int n_unmapped_recovered = 0;
  auto recovery_slot = [&](int vecp_i) {
    return recovery_slots_ok ? sk::LVWORK + vecp_i - first_recovered_row + 1 : 0;
  };

  auto push_particle_words = [&](int vecp_i, int lpa, int lmain,
                                 VecpStatus status, bool chargedPa) {
    // PSHSCT uses the beam spot when no usable primary vertex exists.
    // Preserve the particle, but do not publish that fallback as a valid cut.
    // It also decided whether a charged PA was kept or recovered as neutral.
    const bool missingSelectionReference = sk::IERRBS != 0 &&
        sk::KVTX(16, 1) > 0 && chargedPa;
    std::int32_t word = vecp_i >= 1 ? sk::LVSELE(vecp_i) : -1;
    std::int32_t lock = vecp_i >= 1 ? sk::LVLOCK(vecp_i) : -1;
    if (status == VecpStatus::RecoveredNeutral && sk::IFLODR == 1) {
      const int slot = recovery_slot(vecp_i);
      if (slot > 0) {
        word = sk::LVSELE(slot);
        lock = static_cast<std::int32_t>(
            static_cast<std::uint32_t>(sk::LVLOCK(slot)) |
            (static_cast<std::uint32_t>(word) & 0x80000001u));
      } else {
        word = lock = -1;
        ++n_unmapped_recovered;
      }
    }
    lvlockCol.push_back(missingSelectionReference ? -1 : lock);
    lvseleCol.push_back(word);
    statusCol.push_back(static_cast<std::int32_t>(status));
    chargeCodeCol.push_back(static_cast<std::int32_t>(std::lround(ph::Q(lmain + 8))));
    vecpChargeCol.push_back(vecp_i >= 1 ? sk::VECP(7, vecp_i) : kNaN);
    codeCol.push_back((ph::IQ(lpa + 3) >> 18) & 0x7F);
    idWordCol.push_back(ph::IQ(lpa + 3));
    detCol.push_back(ph::IQ(lpa + 2));
    massCodeCol.push_back(static_cast<std::int32_t>(std::lround(ph::Q(lmain + 10))));
    energyErrorCol.push_back(ph::Q(lmain + 19));
    lengthCol.push_back(ph::Q(lmain + 9));
  };

  // Particle built from VECP row vecp_i: 4-vector, mass and charge as stored.
  auto vecp_particle = [&](int vecp_i, float charge) {
    auto pfo = pfoCol.create();
    pfo.setMomentum({sk::VECP(1, vecp_i), sk::VECP(2, vecp_i), sk::VECP(3, vecp_i)});
    pfo.setEnergy(sk::VECP(4, vecp_i));
    pfo.setMass  (sk::VECP(5, vecp_i));
    pfo.setCharge(charge);
    return pfo;
  };

  // AABTAG's impact parameters, keyed by the PA they belong to. Empty when
  // AABTAG produced nothing usable for this event.
  const auto lpa_to_btag = aabtag::lpaToTrack();

  forEachPA([&](int lpa, int paIdx) {
    // PA.MAIN: per-track summary. Charge code at Q(LMAIN+8):
    //   0 = neutral, 1 = positive, 2 = negative, 3 = undefined.
    const int lmain = lphpa("MAIN", lpa);
    if (lmain <= 0) { mark_pa_skipped(paIdx); return; }
    const int charge_code = static_cast<int>(std::lround(ph::Q(lmain + 8)));

    if (charge_code == 0) {
      // ------- Neutral path -------
      // No PA.TRAC. 4-momentum from VECP. Emit only the Particle.
      const int vecp_i = find_vecp_index(lpa, /*want_charged=*/false);
      if (vecp_i < 1) { mark_pa_skipped(paIdx); return; }

      record_particle(vecp_particle(vecp_i, 0.f), vecp_i, paIdx);
      // d0PV/z0PV/d0BS are parallel to Tracks (charged-only) and so are not
      // pushed for neutrals; the per-particle words are.
      push_particle_words(vecp_i, lpa, lmain, VecpStatus::OwnRow,
                          /*chargedPa=*/false);
      return;
    }

    // ------- Charged path -------
    // PA.TRAC: perigee (d0, z0, theta, phi, 1/R) at +2..+6;
    // 15-element lower-tri weight matrix at +7..+21.
    const int ltrac = lphpa("TRAC", lpa);
    if (ltrac <= 0) { mark_pa_skipped(paIdx); return; }
    const float D0_dp    = ph::Q(ltrac + 2);    // cm
    const float Z0_dp    = ph::Q(ltrac + 3);    // cm
    const float theta_dp = ph::Q(ltrac + 4);    // rad
    const float phi_dp   = ph::Q(ltrac + 5);    // rad
    const float invR_dp  = ph::Q(ltrac + 6);    // 1/cm (signed)
    std::array<float, 15> Wraw{};
    for (int k = 0; k < 15; ++k) Wraw[k] = ph::Q(ltrac + 7 + k);

    const auto helix = Helix::fromPerigee(D0_dp, Z0_dp, theta_dp, phi_dp,
                                          invR_dp, &Wraw);
    if (!helix.valid()) {
      mark_pa_skipped(paIdx);
      return;   // degenerate sin(theta) -> skip
    }

    auto trk = trkCol.create();
    trk.addToTrackStates(helix.toTrackState(edm4hep::TrackState::AtIP));

    // AABTAG measures an impact parameter for the subset of tracks it can
    // use, against its own primary vertex. That is a property of the track,
    // so it rides here as a state at that vertex rather than in a parallel
    // array; a track AABTAG skipped simply has no AtVertex state.
    //
    // D0 is negated into the EDM4hep convention, as the perigee above is
    // (Helix::fromPerigee) -- AABTAG stores the DELPHI sign. Z0 is not
    // negated, matching the same routine. Only these two components are
    // measured; the rest stay NaN rather than zero, which would claim a
    // measurement that was never made.
    if (auto it = lpa_to_btag.find(lpa); it != lpa_to_btag.end()) {
      trk.addToTrackStates(aabtag::vertexState(it->second));
    }

    // Track elements reconstructed from this PA, decoded by
    // TrackElementsWriter. Linked here, while the track is still mutable.
    if (ctx_.track_elements) {
      const auto& te = *ctx_.track_elements;
      if (paIdx < static_cast<int>(te.pa_to_segments.size())) {
        for (const auto& seg : te.pa_to_segments[paIdx]) trk.addToTracks(seg);
      }
      if (paIdx < static_cast<int>(te.pa_to_plane_hits.size())) {
        for (const auto& hit : te.pa_to_plane_hits[paIdx]) {
          trk.addToTrackerHits(hit);
        }
      }
    }

    // Extrapolation states for this PA, decoded by TraxWriter.
    if (ctx_.trax && paIdx < static_cast<int>(ctx_.trax->pa_to_states.size())) {
      for (const auto& st : ctx_.trax->pa_to_states[paIdx]) {
        trk.addToTrackStates(st);
      }
    }

    // chi2 / ndf from PA.MAIN. +26/+27 (with VD) preferred, fallback to
    // +16/+17 (without VD). SKELANA sanitises ndf to [0, 1000]. With
    // neither fit available the track says so: ndf = 0 and chi2 = NaN.
    trk.setChi2(kNotMeasured);
    float chi2_vd = ph::Q(lmain + 26);
    int   ndf_vd  = static_cast<int>(std::lround(ph::Q(lmain + 27)));
    if (ndf_vd < 0 || ndf_vd > 1000) ndf_vd = 0;
    if (ndf_vd > 0 && chi2_vd > 0.f) {
      trk.setChi2(chi2_vd);
      trk.setNdf(ndf_vd);
    } else {
      float chi2_no_vd = ph::Q(lmain + 16);
      int   ndf_no_vd  = static_cast<int>(std::lround(ph::Q(lmain + 17)));
      if (ndf_no_vd < 0 || ndf_no_vd > 1000) ndf_no_vd = 0;
      if (ndf_no_vd > 0) {
        trk.setChi2(chi2_no_vd);
        trk.setNdf(ndf_no_vd);
      }
    }

    // The PA's own charged VECP row: required for the 4-momentum, LVLOCK and
    // QTRAC reads. Without one, the neutral row SKELANA re-entered for this
    // PA (if any) is the calorimeter recovery. No perigee-momentum fallback.
    const int vecp_i = find_vecp_index(lpa, /*want_charged=*/true);
    const int recovered_i = find_vecp_index(lpa, /*want_charged=*/false);

    // Vertex-Detector hits SKELANA assigned to this track, decoded by
    // VdHitsWriter and grouped by the same charged-track ordinal. Linked
    // here, while the track is still mutable.
    if (ctx_.vd_hits && vecp_i >= 1) {
      const auto& vd = ctx_.vd_hits->vecp_to_hits;
      if (vecp_i < static_cast<int>(vd.size())) {
        for (const auto& hit : vd[static_cast<std::size_t>(vecp_i)]) {
          trk.addToTrackerHits(hit);
        }
      }
    }

    // Charge sign from PA.MAIN. Code 3 ("undefined") -> NaN: the track
    // exists but its sign is unknown, which neither +-1 nor 0 can say. The
    // raw code and SKELANA's VECP(7) ride alongside (push_particle_words),
    // so a consumer can choose its own reading.
    float charge = kNotMeasured;
    if      (charge_code == 1) charge = +1.f;
    else if (charge_code == 2) charge = -1.f;

    VecpStatus status = VecpStatus::OwnRow;
    auto pfo = vecp_i >= 1 ? vecp_particle(vecp_i, charge) : pfoCol.create();
    if (vecp_i >= 1) {
      if (recovered_i >= 1) ++n_dual_rows;   // not produced by IFLODR=1
    } else {
      // No charged VECP row: the PA keeps its track, but it has no 4-vector
      // of its own. NaN, never 0, which is a measurable value.
      status = recovered_i >= 1 ? VecpStatus::RecoveredSource
                                : VecpStatus::NoRow;
      ++(recovered_i >= 1 ? n_recovered : n_no_row);
      pfo.setMomentum({kNotMeasured, kNotMeasured, kNotMeasured});
      pfo.setEnergy(kNotMeasured);
      pfo.setMass(kNotMeasured);
      pfo.setCharge(charge);
    }
    pfo.addToTracks(trk);
    record_particle(pfo, vecp_i, paIdx);

    // LVLOCK: per-VECP track-quality bitmask. -1 sentinel if no VECP
    // match. Stored as int32 so bit 32 (REMCLU overlap) is preserved.
    push_particle_words(vecp_i, lpa, lmain, status, /*chargedPa=*/true);

    // BS/PV-corrected impact parameters from sk::QTRAC(38..40, vecp_i), and
    // the perigee words QTRAC(4..23). PSCTRA indexes by the charged-VECP
    // ordinal which matches vecp_i when IFLODR=1 (charged-first ordering).
    push_qtrac_or_nan(vecp_i, ltrac);

    // SKELANA's recovered neutral: its own entry, directly after the PA it
    // came from, and the one place its energy is counted (the source has
    // NaN kinematics). It owns the PA's showers (CalorimeterWriter via
    // calorimeterOwner) and the VECP-row links (truth, STIC, pi0); the
    // source keeps the track and all PA-keyed track information.
    if (status == VecpStatus::RecoveredSource) {
      auto npfo = vecp_particle(recovered_i, 0.f);
      npfo.addToParticles(pfo);
      const int n_idx = static_cast<int>(result.particle_handles.size());
      result.particle_handles.push_back(npfo);
      result.vecp_to_particle[recovered_i] = n_idx;
      if (paIdx >= static_cast<int>(result.pa_to_recovered_neutral.size())) {
        result.pa_to_recovered_neutral.resize(
            static_cast<std::size_t>(paIdx + 1), -1);
      }
      result.pa_to_recovered_neutral[paIdx] = n_idx;
      push_particle_words(recovered_i, lpa, lmain,
                          VecpStatus::RecoveredNeutral, /*chargedPa=*/true);
    }
  });

  if (n_recovered + n_no_row + n_dual_rows > 0) {
    std::cerr << "delphi_edm4hep::tracking: WARNING: charged PFO(s) without a"
              << " charged VECP row in this event: " << n_recovered
              << " recovered by SKELANA as neutral, " << n_no_row
              << " with no row (NaN 4-momentum)";
    if (n_dual_rows > 0) {
      std::cerr << "; " << n_dual_rows << " with both a charged and a"
                << " neutral row (neutral row not emitted)";
    }
    std::cerr << "\n";
  }
  if (n_unmapped_recovered > 0) {
    std::cerr << "delphi_edm4hep::tracking: WARNING: " << n_unmapped_recovered
              << " recovered neutral(s) whose SKELANA work slot could not be"
              << " identified (recovered rows are not the last INREC = "
              << sk::INREC << " VECP rows): SelectionWord/Flag -1\n";
  }

  // Push all collections into the Frame via the base class's put().
  // Handles in `result` remain valid afterwards.
  put(std::move(trkCol),    "TRAC", "Tracks", Provenance::Derived);
  put(std::move(pfoCol),    "MAIN", "Particles", Provenance::Derived);
  put(std::move(lvlockCol), "VECP", "Particles_SelectionFlag", Provenance::Derived);
  put(std::move(lvseleCol), "VECP", "Particles_SelectionWord", Provenance::Derived);
  put(std::move(statusCol), "MAIN", "Particles_VecpStatus",     Provenance::Derived);
  put(std::move(chargeCodeCol), "MAIN", "Particles_ChargeCode", Provenance::Transcribed);
  put(std::move(vecpChargeCol), "VECP", "Particles_Charge",     Provenance::Transcribed);
  put(std::move(codeCol),   "MAIN", "Particles_ReconstructionCode", Provenance::Transcribed);
  put(std::move(idWordCol), "MAIN", "Particles_IdentificationWord", Provenance::Transcribed);
  put(std::move(detCol),    "MAIN", "Particles_DetectorMask",       Provenance::Transcribed);
  put(std::move(massCodeCol),    "MAIN", "Particles_MassCode",      Provenance::Transcribed);
  put(std::move(energyErrorCol), "MAIN", "Particles_EnergyError",   Provenance::Transcribed);
  put(std::move(lengthCol), "MAIN", "Particles_TrackLength",        Provenance::Transcribed);
  put(std::move(d0PvCol),   "QTRAC", "Tracks_d0PV", Provenance::Derived);
  put(std::move(z0PvCol),   "QTRAC", "Tracks_z0PV", Provenance::Derived);
  put(std::move(d0BsCol),   "QTRAC", "Tracks_d0BS", Provenance::Derived);
  for (int w = 0; w < kPerigeeWords; ++w) {
    put(std::move(perigeeCols[w]), "QTRAC", kPerigeeNames[w], Provenance::Derived);
  }
  put(std::move(refittedCol), "QTRAC", "Tracks_PerigeeRefitted", Provenance::Derived);
  put(std::move(fitQualityCol), "QTRAC", "Tracks_FitQuality", Provenance::Derived);
  put(std::move(qtracEnergyErrorCol), "QTRAC", "Tracks_EnergyError", Provenance::Derived);

  // Hand off to downstream writers via the shared context.
  ctx_.tracking = std::move(result);
}

}  // namespace delphi_edm4hep::tracking
