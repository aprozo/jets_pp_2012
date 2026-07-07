# Systematics

Config-as-code systematic band. The variation matrix is a typed, committed list
in `new_ana/config.h` (`Systematic` struct + `Systematics()`); the pipeline runs
each variation and stamps its provenance into every output. No environment flags,
no table file.

## The variation matrix — mirrors the published composition

The sources and sizes follow Dmitry's published band (star-jet
`default.nix:1213-1245`: quadrature of EMC scale, track efficiency, embedding
statistics, UE fraction; luminosity is a separate normalization statement):

| variation | acts in | what it does |
|---|---|---|
| `emcUp/Down` | response | per-jet reco-pT shift ±√(((1−rt)·1.1%)² + (rt·3.2%)²) — BEMC tower scale 3.2%, TPC track scale 1.1%, weighted by the jet's own neutral fraction rt |
| `trkEffUp/Down` | response | 1% track-efficiency equivalent: per-jet ±1%·(1−rt) reco-pT shift (≡ thinning data tracks, opposite sign) |
| `ueUp/Down` | data | detector-side UE-subtraction fraction 1.18 / 0.86 (`pt_corrected + (1−f)·area·ρ`); bypasses the JPX data cache |
| `trigEffUp/Down` | data | ±1.5% on the measured C(pt) — the correction's own measurement precision (plateau fit ±0.8%, turn-on ±1–2%) |
| `unfoldReg` | per-trigger | Bayes nIter 2→3 (regularization dependence) |
| `jpxDamp` | JPX | Tikhonov damping 0→0.030 (unregularized → damped inversion) |
| `embStatUp/Down` | JPX | embedding-statistics: 200 Gaussian resamplings of the coarse response ingredients, canonical = nominal ± 1σ(toys) |

**Luminosity is NOT in the per-bin band** — it is a correlated normalization
uncertainty quoted separately (Dmitry: "10% luminosity uncertainty not shown").
Both analyses use the same Zilong/Dunlop per-run luminosities, so it largely
cancels in me/Dmitry ratios anyway.

`needsResponse()` (jesShift/jerSmear/emcSign/trkSign ≠ 0) decides whether the
response is rebuilt (shape) or the nominal one is reused. The nominal call uses
the `"nominal"` preset — pure identity, the physics pipeline is unchanged.

## Provenance

Every response and xsec ROOT file carries the full variation config as `TNamed`
keys (written by `config.h::StampProvenance`). `combine.C` accepts ONLY files
whose stamped variation is currently in `Systematics()` — solver cross-checks
(`_invert`) and stale files from removed variations never enter the band.

## Run

```bash
bash systematics/run_systematics.sh
```

which (for every trigger in `config.h` AND the JPX combination):
1. builds the nominal response + a rebuilt response per shape variation;
2. runs the cross section for every preset → `xsec_<T>_R0.5_<name>.root` /
   `xsec_JPX_R0.5_<name>.root`;
3. `combine.C` takes the per-bin quadrature envelope →
   `xsec_<T>_R0.5_systband.root` (including `<T>` = JPX; triggers with no
   variation files are skipped).

`list_systematics.C` is the single bridge from `Systematics()` to the shell.
Edit the sources in `config.h`, nowhere else.
