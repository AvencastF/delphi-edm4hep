#ifndef SKELANA_TRACKNTP_HPP
#define SKELANA_TRACKNTP_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
/*  +KEEP,TRACKNTP.  (stdcdes.car)  Variables used in track selection
 *
 *      LVSELE(3*MTRACK) - bitted word per track (old LVLOCK + track
 *                         selection cuts + REMCLU related locking),
 *                         COMMON /TRKINF/, indexed like LVLOCK by VECP row
 *
 *  Bit meanings, PSHORTSEL header (skelana.car), DELPHI numbering (bit 1 =
 *  least significant):
 *      1      standard selection (global)
 *      2-11   cut that caused the rejection: ICTFAIL (charged, PSHSCT) or
 *             INTFAIL (neutral, PSHSNT) bits 1-10 shifted up by one
 *      20-22  recovery step result (20 calorimeter, 21 refit with beam-spot
 *             constraint, 22 mammoth)
 *      30     REMCLU cluster added to VECP
 *      31     PA in REMCLU cluster
 *      32     multivertex structure handling
 *
 *  The other TRACKNTP commons (/CHGEDTRK/ ... ICTFAIL, /NEUTRK/ ... INTFAIL,
 *  /RECOVER/ INREC) hold scalars of the track last examined, not per-track
 *  arrays, and are not bound here.
 *
 *  Not vendored from delphi-nanoaod: added by delphi-edm4hep (2026-11-12).
 */
    extern "C" struct
    {
        int lvsele[3 * MTRACK];
    } trkinf_;

    inline int &LVSELE(int i) { return trkinf_.lvsele[i - 1]; }
} // namespace skelana

#endif // SKELANA_TRACKNTP_HPP
