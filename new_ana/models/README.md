# models: generators against the data

Particle-level inclusive jet spectra of Pythia 6 (Perugia 2012, the embedding generator), Pythia 8
(Monash and the RHIC Detroit tune) and Herwig 7.3 (default tune), with the Stage-1 particle-level
jet definition, compared with the unfolded JP1 cross section; and the hadronisation correction
C_had = sigma_particle / sigma_parton for the fixed-order comparison. Runs on the host in the LCG
view (`../../site.sh`).

```bash
bash new_ana/models/run_models.sh gen pythia8 100000   # one generator, 13 pT-hat bins
bash new_ana/models/run_models.sh chad                 # Perugia 2012 and its five variations, particle and parton level
bash new_ana/models/run_models.sh spectra              # jets from the trees -> results/models/spectra_<generator>.root
bash new_ana/models/run_models.sh compare              # tables against the data -> results/models/models_R<R>_{ue,noue}.{txt,root}, chad_R0.5.{txt,root}
bash new_ana/models/run_models.sh all
```

The five small programs are compiled into `bin/` (not tracked) by the driver itself, whenever a
source or `tree.h` is newer than its binary; `run_models.sh build` forces a rebuild.

## Changing a parameter

Every physics setting of a generator is written out, with a comment per line, at the top of its
program, so the whole configuration can be read off the source:

| file | generator | where the settings are |
|---|---|---|
| `gen_pythia6.cc` | Pythia 6.4, Perugia 2012, the tune of the embedding request | the block marked `settings`: tune, PARP(90), process, energy, the 13 undecayed species |
| `gen_pythia8_monash.cc` | Pythia 8, default Monash 2013 | the list `kSettings`, in Pythia 8's own `Key = value` language |
| `gen_pythia8_detroit.cc` | Pythia 8, Detroit tune, PRD 105 (2022) 016011 Table III | the same, with the tune block |
| `herwig.in` | Herwig 7.3, default tune | Herwig's own input file, `@PTMIN@`, `@PTMAX@`, `@OUT@` substituted per bin |

To change a setting, a tune parameter, the beam energy, which species stay undecayed, edit it there
and rerun `gen` for that generator; the driver recompiles a program whose source is newer than its
binary. What is not in the sources is what changes from job to job: the pT-hat bin, the number of
events, the random seed and, for the C_had variations, the PYTUNE index, which `run_models.sh` passes
on the command line. Every tree records its settings in a `settings` TNamed:

```bash
root -l ../../output/models/pythia8/pt11_15.root -e 'cout << ((TNamed*)gFile->Get("settings"))->GetTitle();'
```

## Steps

1. `gen_pythia6.cc`, `gen_pythia8_monash.cc`, `gen_pythia8_detroit.cc`, `herwig.in` + `hepmc2tree.cc` — one particle tree per
   pT-hat bin (`tree.h`: final state with |eta| < 3, the hard-process pT, the bin cross section), the
   13 bins of the embedding request, the same 13 particles left undecayed. Trees in
   `../../output/models/`.
2. `jetspec.cc` — anti-kT jets with the off-axis-cone UE subtraction at R = 0.2-0.5, weighted by
   sigma / N per bin, the outlier filter of the published analysis per bin; Pythia 6 also with the soft
   reweight and the soft-sample corrections of the embedding.
3. `compare.C` — the spectra rebinned to the data bins against JP1 (with its band), the combined
   levels, the min-bias level and, at R = 0.5, Table III. Generators are quoted from 9.7 GeV
   (10 GeV on the 5-60 GeV bins).
4. `chad.C` — C_had per bin from the Perugia 2012 samples at particle and parton level (same
   events); uncertainty from the radiation (371/372), fragmentation (376/377) and Innsbruck (373)
   variations in quadrature.

Figures: `bash new_ana/plots.sh` (`analysis/plot.C`).
