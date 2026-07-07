# Systematics

Config-as-code systematic band. The variation matrix is a typed, committed list
in `new_ana/config.h` (`Systematic` struct + `Systematics()`); the pipeline runs
each variation and stamps its provenance into every output. No environment flags,
no table file.

## How a variation is defined

`config.h::Systematics()` returns a list of `Systematic` presets:

```cpp
inline std::vector<Systematic> Systematics() {
   return {
      {"nominal"},
      SystJES("jesUp", +0.03), SystJES("jesDown", -0.03),   // reco energy scale (shape)
      SystJER("jerUp", 0.05),                               // reco resolution   (shape)
      SystNorm("trigEffUp", 1.0, 1.05), SystNorm("trigEffDown", 1.0, 0.95),
      SystNorm("lumiUp", 1.086, 1.0),   SystNorm("lumiDown", 0.914, 1.0),
      SystIter("unfoldReg", 3),                             // Bayes nIter 2->3
   };
}
```

Each field shifts exactly one thing:

| field          | acts in            | kind                                    |
|----------------|--------------------|-----------------------------------------|
| `jesShift`     | `unfold.cxx`       | SHAPE — response reco energy scale; response is rebuilt |
| `jerSmear`     | `unfold.cxx`       | SHAPE — extra Gaussian reco smear; response is rebuilt  |
| `trigEffScale` | `cross_section.cpp`| normalization — flat scale on C(pt)     |
| `lumiScale`    | `cross_section.cpp`| normalization — luminosity scale        |
| `nIter`        | `cross_section.cpp`| unfolding — Bayes iterations (reuses nominal response) |

`needsResponse()` (jesShift/jerSmear ≠ 0) decides whether the response is rebuilt
(shape) or the nominal one is reused (everything else). The nominal call
`unfold()` / `cross_section()` uses the `"nominal"` preset — pure identity, so the
physics pipeline is unchanged.

## Provenance

Every response and xsec ROOT file carries the full variation config as `TNamed`
keys (`variation`, `jesShift`, `jerSmear`, `lumiScale`, `trigEffScale`, `nIter`),
written by `config.h::StampProvenance`. "What produced this histogram" is
answerable from the file itself.

## Run

```bash
bash systematics/run_systematics.sh
```

which (for every trigger in `config.h`):
1. builds the nominal response + a rebuilt response per shape variation;
2. runs `cross_section` for every preset → `xsec_<T>_R0.5_<name>.root`;
3. `combine.C` takes the per-bin up/down envelope → `xsec_<T>_R0.5_systband.root`
   (`canonical` / `systematic` / `reference`).

`list_systematics.C` is the single bridge from `Systematics()` to the shell
(prints the names); edit the sources in `config.h`, nowhere else.

## Notes

- The JES/JER values here are illustrative — set them to the measured 1σ
  uncertainties.
- Luminosity is a correlated NORMALIZATION uncertainty; many analyses quote it
  separately rather than folding it into the per-bin band. Drop the `lumi*` rows
  from `Systematics()` if you prefer that convention.
- The unfolding systematic can also be taken as the Bayes vs matrix-inversion
  difference — run `cross_section_inverse` and add its `xsec_<T>_R0.5_invert.root`
  as a variation (rename it to `xsec_<T>_R0.5_<name>.root` so `combine.C` picks
  it up).
