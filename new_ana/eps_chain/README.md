# eps_chain: the min-bias chain efficiency

The min-bias level is normalised by the sampled VPDMB-nobsmd luminosity times
eps_chain = P(VPDMB fired and |vz_vpd - vz| < 6 cm | a jet event), measured in the jet-patch data:
jets of the JP0, JP1 and JP2 samples, the min-bias hardware bit and the VPD match read from the pico
header. The result is radius-independent within errors; `stage2/unfold.C` uses the error-weighted
JP0 constant.

Rerun only after a new Stage-1 data production; the current result is tracked in `../inputs/`.

```bash
bash new_ana/eps_chain/run_bits.sh [nchunks] [njobs]         # pass 1 over the picos -> bits4/
# inside the container, from this directory, per radius:
root -l -b -q 'vpd_trigger_eff_jp_pass2d.C("bits4","0.5")'   # pass 2 -> bits4/eps_chain3_R0.5.root
```

1. `vpd_trigger_eff_jp_pass1b.C` — per jet-patch-fired event: the hardware min-bias ids, the JP ids,
   the VPD flags and the join key of the Stage-1 tree.
2. `vpd_trigger_eff_jp_pass2d.C` — joins pass 1 to the Stage-1 tree, fills denominator and chain
   histograms per trigger sample and jet pT (prescale-weighted), writes `eps_chain3_R<R>.root`;
   `run_bits.sh` copies it to `../inputs/`.
