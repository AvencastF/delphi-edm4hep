#!/usr/bin/env python3
"""Alignment / relation audit for delphi_edm4hep EDM4hep output.

A companion UserData array is index-aligned with an object collection: row i
describes element i of that collection. This scanner checks that contract on a
real file, plus the RecDqdx -> Track relation, the primary-vertex flags and the
MAIN VecpStatus contract (NaN exactly where no own VECP row exists, each
recovered neutral directly after its one source, every shower owned by at most
one particle, charge consistent with the raw charge code, VECP charge NaN
exactly where no own VECP row exists, ReconstructionCode = bits 19-25 of the
PA identification word, recovered neutrals repeating their source's PA words,
SKELANA perigee words NaN exactly where no own VECP row exists and equal to
the AtIP state wherever they are not marked refitted, the QTRAC fit words NaN
exactly there too, an unrefitted QTRAC energy error equal to PA.MAIN word 19
and a refit's fit quality a probability >= 0.01, the LVSELE selection word -1
exactly without an own row, SelectionFlag equal to its bits 1 and 32, and bit
20 on every recovered neutral). It
guards the bug class where a per-particle array is emitted alongside a track
collection, so consumers read the wrong row. Value-level checks cannot see it,
and it only shows on events with a reconstructed neutral (Particles > Tracks):
a clean two-track event has Particles == Tracks and hides it.

Usage:
    align_audit.py <file.edm4hep.root>
    DELPHI_EDM4HEP_SAMPLE=<file.edm4hep.root> align_audit.py

Exit 0 if all contracts hold, 77 if no input is configured (CTest skip), and
1 if an explicitly configured input is missing or any contract fails.
"""
import math
import os
import struct
import sys

# A companion array is index-aligned with an object collection. The array is
# either named after its parent plus a suffix, or carries the parent's kind as
# a token.
PARENT_BY_KIND = {
    "Tracks":    "TRAC_Tracks",
    "Particles": "MAIN_Particles",
}


def parent_of(name, names):
    """The collection a companion array runs parallel to, or None."""
    stem = name.rsplit("_", 1)[0]
    if stem in names:
        return stem
    tokens = name.split("_")
    for token in tokens[1:]:
        if token in PARENT_BY_KIND:
            parent = f"{tokens[0]}_{PARENT_BY_KIND[token]}"
            return parent if parent in names else None
    return None


def same(a, b):
    """Equal, NaN equal to NaN."""
    return a == b or (a != a and b != b)


def f32(x):
    """x rounded to float32, as the converter computes it."""
    return struct.unpack("f", struct.pack("f", x))[0]


def perigee_violations(perigee, track, status, fit=None, pa_error=None):
    """SKELANA perigee and fit words of one track against its status, its
    particle's PA energy error and its AtIP state."""
    k = track.getObjectID().index
    flag = int(perigee["Refitted"][k])
    words = [float(v[k]) for n, v in perigee.items() if n != "Refitted"]
    bad = []
    if flag not in (-1, 0, 1) or (flag == -1) != (status in (1, 3)):
        bad.append(f"PerigeeRefitted {flag} with status {status}")
    if any(w != w for w in words) != (flag == -1) or \
            all(w != w for w in words) != (flag == -1):
        bad.append(f"perigee words NaN pattern wrong for PerigeeRefitted {flag}")
    if fit:
        # QTRAC(28)/(41): NaN exactly without an own row; unrefitted QTRAC(41)
        # is PA.MAIN word 19; a refit's QTRAC(28) is CONFPV's chi2
        # probability, which CONFPV only accepts at >= 0.01.
        quality, error = float(fit["FitQuality"][k]), float(fit["EnergyError"][k])
        if (quality != quality) != (flag == -1) or (error != error) != (flag == -1):
            bad.append(f"fit words NaN pattern wrong for PerigeeRefitted {flag}")
        if flag == 0 and pa_error is not None and not same(error, pa_error):
            bad.append("unrefitted QTRAC EnergyError differs from PA.MAIN word 19")
        if flag == 1 and not f32(0.01) <= quality <= 1:
            bad.append("refitted FitQuality outside the CONFPV probability range")
    if flag == 0:
        # Not refitted: the words are PA.TRAC, the source of the AtIP state
        # (Helix::fromPerigee: D0 = -10 d0, Z0 = 10 z0, omega = 1/R / 10).
        st = list(track.getTrackStates())[0]
        g = {n: float(perigee[n][k]) for n in ("D0", "Z0", "Theta", "Phi", "InvR")}
        want = (f32(-10 * g["D0"]), f32(10 * g["Z0"]), g["Phi"],
                f32(g["InvR"] / 10), 1 / math.tan(g["Theta"]))
        got = (st.D0, st.Z0, st.phi, st.omega, st.tanLambda)
        tol = (1e-9, 1e-9, 1e-9, 1e-12, 1e-5)  # tanLambda: float cos/sin
        if int(st.location) != 1 or not all(
                math.isclose(a, b, rel_tol=1e-5, abs_tol=t)
                for a, b, t in zip(want, got, tol)):
            bad.append("unrefitted perigee words differ from the AtIP state")
    return bad


def vecp_status_violations(fr, names, prefix="sDST"):
    """VecpStatus contract violations in one frame, as a list of strings."""
    main, stat = f"{prefix}_MAIN_Particles", f"{prefix}_MAIN_Particles_VecpStatus"
    if main not in names or stat not in names:
        return []
    ps, st = fr.get(main), fr.get(stat)
    if len(ps) != len(st):
        return []  # reported by the length check
    # Charge words, checked when present (absent before they were added).
    code, vecq = f"{main}_ChargeCode", f"{prefix}_VECP_Particles_Charge"
    cc = vc = None
    if code in names and vecq in names:
        cc, vc = fr.get(code), fr.get(vecq)
        if len(cc) != len(ps) or len(vc) != len(ps):
            cc = vc = None  # reported by the length check
    # PA words and SKELANA perigee words, checked when present (absent before
    # they were added) and when their lengths match (else the length check
    # reports them).
    pa = {w: fr.get(f"{main}_{w}") for w in
          ("IdentificationWord", "MassCode", "EnergyError",
           "ReconstructionCode") if f"{main}_{w}" in names}
    if any(len(v) != len(ps) for v in pa.values()):
        pa = {}
    tracks = fr.get(f"{prefix}_TRAC_Tracks") \
        if f"{prefix}_TRAC_Tracks" in names else []
    stem = f"{prefix}_QTRAC_Tracks_Perigee"
    perigee = {n[len(stem):]: fr.get(n) for n in names if n.startswith(stem)}
    if "Refitted" not in perigee or \
            any(len(v) != len(tracks) for v in perigee.values()):
        perigee = {}
    fit = {w: fr.get(f"{prefix}_QTRAC_Tracks_{w}") for w in
           ("FitQuality", "EnergyError")
           if f"{prefix}_QTRAC_Tracks_{w}" in names}
    if len(fit) != 2 or any(len(v) != len(tracks) for v in fit.values()):
        fit = None
    # SKELANA's LVSELE word and the LVLOCK verdict built from its bits 1, 32.
    sel = {w: fr.get(f"{prefix}_VECP_Particles_{w}") for w in
           ("SelectionWord", "SelectionFlag")
           if f"{prefix}_VECP_Particles_{w}" in names}
    if len(sel) != 2 or any(len(v) != len(ps) for v in sel.values()):
        sel = None
    bad = []
    owners = {}
    for i, p in enumerate(ps):
        s = int(st[i])
        v = p.getMomentum()
        finite = all(x == x for x in (v.x, v.y, v.z, p.getEnergy(), p.getMass()))
        if s not in (0, 1, 2, 3):
            bad.append(f"status {s} undefined")
        if finite == (s in (1, 3)):
            bad.append(f"status {s} with finite={finite}")
        if s in (1, 3) and len(p.getTracks()) != 1:
            bad.append(f"status {s} without its track")
        if s == 1 and (i + 1 >= len(ps) or int(st[i + 1]) != 2):
            bad.append("recovered source not followed by its neutral")
        if s == 2:
            src = [q.getObjectID().index for q in p.getParticles()]
            if src != [i - 1] or int(st[i - 1]) != 1:
                bad.append("recovered neutral not linked to the preceding source")
            if p.getCharge() != 0 or len(p.getTracks()):
                bad.append("recovered neutral carries charge or a track")
        if cc is not None:
            code, vq, q = int(cc[i]), float(vc[i]), p.getCharge()
            # charge derives from the raw code; a recovered neutral repeats
            # its source's code but is neutral by construction.
            want = 0.0 if s == 2 else {0: 0.0, 1: 1.0, 2: -1.0}.get(code)
            if code not in (0, 1, 2, 3):
                bad.append(f"charge code {code} undefined")
            elif not (q == want or (want is None and q != q)):
                bad.append(f"charge {q} inconsistent with charge code {code}")
            if (vq != vq) != (s in (1, 3)):
                bad.append(f"VECP charge {vq} with status {s}")
            if s == 2 and (vq != 0 or code == 0 or code != int(cc[i - 1])):
                bad.append("recovered neutral VECP charge or code wrong")
        if "IdentificationWord" in pa and "ReconstructionCode" in pa and \
                (int(pa["IdentificationWord"][i]) >> 18 & 0x7F) != \
                int(pa["ReconstructionCode"][i]):
            bad.append("ReconstructionCode is not IdentificationWord bits 19-25")
        if s == 2 and any(not same(v[i], v[i - 1]) for v in pa.values()):
            bad.append("recovered neutral PA words differ from its source")
        if sel is not None:
            word = int(sel["SelectionWord"][i]) & 0xFFFFFFFF
            flag = int(sel["SelectionFlag"][i]) & 0xFFFFFFFF
            if (word == 0xFFFFFFFF) != (s in (1, 3)):
                bad.append(f"SelectionWord -1 pattern wrong for status {s}")
            elif s in (0, 2) and flag != 0xFFFFFFFF and \
                    flag != word & 0x80000001:
                bad.append("SelectionFlag is not SelectionWord bits 1 and 32")
            if s == 2 and not word >> 19 & 1:
                bad.append("recovered neutral without SelectionWord bit 20")
        if perigee and len(p.getTracks()) == 1:
            bad.extend(perigee_violations(
                perigee, list(p.getTracks())[0], s, fit,
                float(pa["EnergyError"][i]) if "EnergyError" in pa else None))
        for c in p.getClusters():
            oid = c.getObjectID()
            key = (oid.collectionID, oid.index)
            if key in owners:
                bad.append("shower attached to two particles")
            owners[key] = i
    return bad


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] \
        else os.environ.get("DELPHI_EDM4HEP_SAMPLE", "")
    if not fn:
        print("align_audit: no input (pass a file or set DELPHI_EDM4HEP_SAMPLE)"
              " -> skipping")
        return 77
    if not os.path.exists(fn):
        print(f"align_audit: configured input not found: {fn}")
        return 1

    from podio import root_io
    reader = root_io.Reader(fn)

    len_fail = {}       # (array, parent) -> [(event, got, expected), ...]
    audited = set()
    unresolved = set()
    dq_missing = dq_total = nev = 0
    dummy_first = dummy_published_as_primary = dummy_chain_marked_primary = 0
    status_fail = {}    # violation -> [event, ...]
    for i, fr in enumerate(reader.get("events")):
        nev += 1
        names = set(fr.getAvailableCollections())

        for name in names:
            coll = fr.get(name)
            if "UserData" not in str(coll.getTypeName()):
                continue
            parent = parent_of(name, names)
            # An array whose parent cannot be resolved is a failure, not a
            # skip: a renamed collection would otherwise disable its own check.
            if parent is None:
                unresolved.add(name)
                continue
            audited.add(name)
            got, expected = coll.size(), fr.get(parent).size()
            if got != expected:
                len_fail.setdefault((name, parent), []).append(
                    (i, got, expected))

        for n in names:
            if n.endswith("_RecDqdx"):
                for dq in fr.get(n):
                    dq_total += 1
                    if dq.getTrack().getObjectID().index < 0:
                        dq_missing += 1

        for v in vecp_status_violations(fr, names):
            status_fail.setdefault(v, []).append(i)

        status_name = "sDST_PV_Vertices_StatusBits"
        primary_name = "sDST_PV_PrimaryVertex"
        if status_name in names and primary_name in names:
            statuses = fr.get(status_name)
            if statuses.size() and (int(statuses[0]) & 0x01):
                dummy_first += 1
                if fr.get(primary_name).size():
                    dummy_published_as_primary += 1
                vertices = fr.get("sDST_PV_Vertices")
                if vertices.size() and vertices[0].isPrimary():
                    dummy_chain_marked_primary += 1

    print(f"align_audit: {fn}  ({nev} events)")
    ok = True
    if nev == 0:
        ok = False
        print("  FAIL  input contains no event frames")
    if unresolved:
        ok = False
        for name in sorted(unresolved):
            print(f"  FAIL  {name} has no parent collection to align against")
    if len_fail:
        ok = False
        for (array, parent), lst in sorted(len_fail.items()):
            e = lst[0]
            print(f"  FAIL  len({array}) != len({parent}) in {len(lst)}/{nev} "
                  f"events (e.g. evt{e[0]}: {e[1]} vs {e[2]})")
    elif audited:
        print(f"  OK    all {len(audited)} companion arrays match their parent")
    if dq_missing:
        ok = False
        print(f"  FAIL  {dq_missing}/{dq_total} RecDqdx rows missing a Track link")
    elif dq_total:
        print(f"  OK    all {dq_total} RecDqdx rows carry a Track link")
    if status_fail:
        ok = False
        for v, evs in sorted(status_fail.items()):
            print(f"  FAIL  VecpStatus: {v} in {len(evs)}/{nev} events "
                  f"(e.g. evt{evs[0]})")
    else:
        print("  OK    MAIN VecpStatus contract holds")
    if dummy_published_as_primary:
        ok = False
        print(f"  FAIL  {dummy_published_as_primary}/{dummy_first} events with "
              "a dummy first DELPHI vertex published sDST_PV_PrimaryVertex")
    elif dummy_first:
        print(f"  OK    all {dummy_first} dummy first DELPHI vertices are "
              "excluded from sDST_PV_PrimaryVertex")
    if dummy_chain_marked_primary:
        ok = False
        print(f"  FAIL  {dummy_chain_marked_primary}/{dummy_first} dummy first "
              "DELPHI vertices remain marked primary in sDST_PV_Vertices")
    elif dummy_first:
        print(f"  OK    all {dummy_first} dummy first DELPHI vertices are marked "
              "non-primary in sDST_PV_Vertices")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
