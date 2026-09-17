# inputs: the numbers Stage-2 needs besides the trees

| file | content | used by |
|---|---|---|
| `lumi_zilong_full.root` | per-run sampled luminosity of every jet-patch and high-tower trigger (`luminosity_JP0`, `_JP1`, `_JP2`, `_HT2`, bin labels = run) | `stage2/unfold.C`, run selection of `build_data.C` / `build_resp.C` |
| `lumi_VPDMB_true.root` | per-run min-bias luminosity (`luminosity_MBtrue`): ZDC luminosity times live time over prescale, effective cross section 3.28 mb | the min-bias level |
| `eps_chain3_R<R>.root` | the chain efficiency per trigger sample and jet pT (`den<s>`, `chain<s>`, `raw<s>`, s = JP0, JP1, JP2), from `../eps_chain/` | `stage2/unfold.C` |
| `jet_cross_section_publishedR0.5.root` | the published Table III (`crossSection_statistic`, `crossSection_systematic`) | every comparison and figure |

The luminosity basis is the per-run sampled luminosity for every trigger; the min-bias luminosity is
never N_events / 25 mb.
