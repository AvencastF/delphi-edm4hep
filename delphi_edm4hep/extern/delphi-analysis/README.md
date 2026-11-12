# extern/delphi-analysis — vendored PHDST / SKELANA headers

Header-only C++ wrappers around the DELPHI PHDST and SKELANA Fortran
common blocks, vendored from the `delphi-analysis` package of
[`delphi-nanoaod`](https://github.com/DickyChant/delphi-nanoaod).

- **Upstream repo:** https://github.com/DickyChant/delphi-nanoaod
- **Vendored from commit:** `fde7b95f6ce8533f58825685980add34dd5e0faf`
  (branch `feature/sim-truth-pv`, vendored 2026-08-19)
- **Upstream path:** `delphi-analysis/include/{phdst,skelana}/`

This replaces the former `extern/delphi-nanoaod` git submodule. Only the
headers this converter actually includes (plus their transitive includes)
are kept — `pscbtg.hpp` was added when b-tagging output was implemented — the umbrella headers (`phdst.hpp`, `skelana.hpp`), the SKELANA
common blocks we never read, and everything else in delphi-nanoaod
(nanoaod writers, python, scripts) are intentionally not vendored.

`skelana/trackntp.hpp` (the `TRKINF` common, SKELANA's per-track selection
word `LVSELE`, and the `RECOVER` counter `INREC`, from the `TRACKNTP` keep in
`stdcdes.car`, plus `LVWORK` from keep `PSCVECJJ`) has no upstream
counterpart: it was written here, in the same style, and a re-copy must keep it.

## Updating

Re-copy the needed headers from a delphi-nanoaod checkout and update the
commit hash above:

```sh
cp <nanoaod>/delphi-analysis/include/phdst/{functions,phciii,uxcom,uxlink}.hpp   phdst/
cp <nanoaod>/delphi-analysis/include/skelana/<needed>.hpp                        skelana/
```

If you add an `#include "phdst/..."` or `#include "skelana/..."` to the
converter that is not vendored here, copy that header (and anything it
includes, e.g. `skelana/mtrack.hpp`) too — the build fails loudly otherwise.

## License

delphi-nanoaod carries no explicit license (see the top-level README's
License section); these files are redistributed here with the upstream
author's knowledge, under the same "all rights reserved until clarified"
status as the rest of this repository.
