#!/usr/bin/env python3
"""Alignment / relation audit for delphi_edm4hep EDM4hep output.

A companion UserData array is index-aligned with an object collection: row i
describes element i of that collection. This scanner checks that contract on a
real file, plus the RecDqdx -> Track relation and the primary-vertex flags. It
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
import os
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
