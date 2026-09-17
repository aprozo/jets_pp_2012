#!/bin/bash
# run_models.sh — the generator comparison of new_ana/models, from the host in the LCG view
# (Pythia 6/8, Herwig 7, ROOT, fastjet, LHAPDF).
# Every physics setting of a generator is written out at the top of its program:
#   gen_pythia6.cc          Pythia 6.4, the tune and settings of the STAR embedding request (Perugia 2012)
#   gen_pythia8_monash.cc   Pythia 8, default Monash 2013 tune
#   gen_pythia8_detroit.cc  Pythia 8, RHIC Detroit tune
#   herwig.in               Herwig 7, default tune (Herwig's own input file)
# To change a parameter, edit it there and rerun gen; the programs in bin/ are recompiled by the driver
# whenever a source or tree.h is newer than the binary.
#   bash new_ana/models/run_models.sh gen <generator> [N]      the 13 pT-hat bins of the embedding request, N events each
#                                                              (default 100000, herwig 20000), 4 bins at a time at nice 19;
#                                                              generators: pythia6 pythia8 pythia8detroit herwig
#   bash new_ana/models/run_models.sh chad [N]                 the hadronisation correction samples: Pythia6 Perugia 2012 and its
#                                                              variants 371 372 373 376 377, particle and parton level of the same
#                                                              events (N default 30000) -> output/models/pythia6_t<tune>[_partons]/
#   bash new_ana/models/run_models.sh spectra [generator ...]  jetspec -> new_ana/results/models/spectra_<generator>.root
#                                                              (the chad samples at R = 0.5 only)
#   bash new_ana/models/run_models.sh compare [R ...]          compare.C -> new_ana/results/models/models_R<R>_{ue,noue}.{txt,root};
#                                                              chad.C -> new_ana/results/models/chad_R0.5.{txt,root}
#                                                              (the figures are drawn by new_ana/plots.sh)
#   bash new_ana/models/run_models.sh build                    force a rebuild of the five programs -> new_ana/models/bin/
#   bash new_ana/models/run_models.sh all                      gen x 4, chad, spectra, compare (radii 0.2 0.3 0.4 0.5)
# Trees: output/models/<generator>/pt<lo>_<hi>.root (gitignored), one log per bin next to them.
set -e

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
source "$repo/site.sh"
source $LCG/setup.sh

# the four generators of the model comparison, and the Perugia 2012 tune variants of C_had
GENS="pythia6 pythia8 pythia8detroit herwig"
CHAD_TUNES="370 371 372 373 376 377"

# the pT-hat bins of the 2023 embedding request, in GeV; -1 = open upper edge
BINS="2 3
3 4
4 5
5 7
7 9
9 11
11 15
15 20
20 25
25 35
35 45
45 55
55 -1"

# bins generated at the same time
NPAR=4

# Compile one of the five executables against the LCG view.
compile() {
   mkdir -p "$here/bin"
   case $1 in
      gen_pythia6)
         g++ -O2 -o "$here/bin/gen_pythia6" "$here/gen_pythia6.cc" $(root-config --cflags --libs) \
             -L$LCG/lib -lpythia6 -lpythia6_dummy -L$(lhapdf-config --libdir) -lLHAPDF -lgfortran
         ;;
      gen_pythia8_monash|gen_pythia8_detroit)
         g++ -O2 -o "$here/bin/$1" "$here/$1.cc" $(root-config --cflags --libs) \
             -I$(pythia8-config --includedir) $(pythia8-config --libs)
         ;;
      hepmc2tree)
         g++ -O2 -o "$here/bin/hepmc2tree" "$here/hepmc2tree.cc" $(root-config --cflags --libs) \
             $(HepMC3-config --cxxflags --libs)
         ;;
      jetspec)
         g++ -O2 -o "$here/bin/jetspec" "$here/jetspec.cc" $(root-config --cflags --libs) \
             $(fastjet-config --cxxflags --libs)
         ;;
   esac
   echo "built bin/$1"
}

# Compile the programs a step needs, if the binary is missing or older than its source or tree.h.
# Called before every step that runs one, so the build never has to be asked for.
need() {
   local prog
   for prog in "$@"; do
      if [ ! -x "$here/bin/$prog" ] || [ "$here/$prog.cc" -nt "$here/bin/$prog" ] ||
         [ "$here/tree.h" -nt "$here/bin/$prog" ]; then
         compile "$prog"
      fi
   done
}

# Force a rebuild of all five, whatever the timestamps say.
build() {
   compile gen_pythia6
   compile gen_pythia8_monash
   compile gen_pythia8_detroit
   compile hepmc2tree
   compile jetspec
}

# The program that makes the trees of one generator.
generator_program() {
   case $1 in
      herwig)         echo hepmc2tree ;;
      pythia8)        echo gen_pythia8_monash ;;
      pythia8detroit) echo gen_pythia8_detroit ;;
      *)              echo gen_pythia6 ;;
   esac
}

# One pT-hat bin of one generator: bin <generator> <N> <lo> <hi> <index>.
# Called through "bash $0 bin ..." by gen(), so it is also a case below.
bin() {
   local g=$1 n=$2 lo=$3 hi=$4 idx=$5
   local tag="pt${lo}_${hi}"
   local out="$repo/output/models/$g"
   # one seed series per generator; all C_had tunes share theirs, so they generate the same events
   local seed
   case $g in
      pythia6)        seed=$((1000 + idx)) ;;
      pythia8)        seed=$((2000 + idx)) ;;
      herwig)         seed=$((3000 + idx)) ;;
      pythia8detroit) seed=$((4000 + idx)) ;;
      pythia6_t*)     seed=$((10000 + idx)) ;;
      *)
         echo "unknown generator $g"
         exit 1
         ;;
   esac
   mkdir -p "$out"
   local log="$out/$tag.log"
   case $g in
      herwig)
         # Herwig writes HepMC in its own run directory; hepmc2tree turns it into the particle tree
         local d="$out/$tag.work"
         local hi2=$hi
         [ "$hi" = "-1" ] && hi2=1000000
         rm -rf "$d"
         mkdir -p "$d"
         sed -e "s/@PTMIN@/$lo/" -e "s/@PTMAX@/$hi2/" -e "s/@OUT@/herwig.hepmc/" "$here/herwig.in" > "$d/herwig.in"
         (cd "$d" && Herwig read herwig.in --repo=$LCG/share/Herwig/HerwigDefaults.rpo -I $LCG/share/Herwig && nice -n 19 Herwig run herwig.run -N $n -s $seed) > "$log" 2>&1
         "$here/bin/hepmc2tree" "$d/herwig.hepmc" "$out/$tag.root" "Herwig 7.3.0 default tune, MEQCD2to2, kT $lo-$hi GeV, pp 200 GeV, 13 particles undecayed, seed $seed" >> "$log" 2>&1
         cat "$d"/herwig-S*.out >> "$log"
         rm -rf "$d"
         ;;
      pythia8|pythia8detroit)
         nice -n 19 "$here/bin/$(generator_program $g)" $lo $hi $n $seed "$out/$tag.root" > "$log" 2>&1
         ;;
      pythia6)
         nice -n 19 "$here/bin/gen_pythia6" $lo $hi $n $seed "$out/$tag.root" > "$log" 2>&1
         ;;
      pythia6_t*)
         # a C_had tune: the same events are written twice, partons before hadronisation and particles after
         mkdir -p "${out}_partons"
         nice -n 19 "$here/bin/gen_pythia6" $lo $hi $n $seed "$out/$tag.root" "${g#pythia6_t}" \
              "${out}_partons/$tag.root" > "$log" 2>&1
         ;;
   esac
   echo "$g $tag done"
}

# Every pT-hat bin of one generator, NPAR bins at a time: gen <generator> [N]
gen() {
   local g=$1 n=$2
   [ -z "$n" ] && { n=100000; [ "$g" = herwig ] && n=20000; }
   need "$(generator_program $g)"
   echo "$BINS" | awk '{print $1, $2, NR}' | xargs -P $NPAR -L 1 bash "$0" bin "$g" "$n"
   ls "$repo/output/models/$g"/pt*.root | wc -l | xargs echo "$g: trees"
}

# The hadronisation-correction samples: the Perugia 2012 tune and its five variants.
chad() {
   local n=${1:-30000}
   need gen_pythia6
   for t in $CHAD_TUNES; do
      gen pythia6_t$t $n
   done
}

# Particle-level jet spectra of the generators that have trees: spectra [generator ...]
spectra() {
   need jetspec
   mkdir -p "$repo/new_ana/results/models"
   local list=${@:-$GENS $(for t in $CHAD_TUNES; do echo pythia6_t$t pythia6_t${t}_partons; done)}
   for g in $list; do
      local radii="0.2 0.3 0.4 0.5"
      # the C_had samples are only needed at the published radius
      [ "${g#pythia6_t}" != "$g" ] && radii="0.5"
      [ -d "$repo/output/models/$g" ] || { echo "spectra: no trees for $g"; continue; }
      "$here/bin/jetspec" "$g" "$repo/output/models/$g" "$repo/new_ana/results/models/spectra_$g.root" "$radii"
   done
}

# The generators against the unfolded data, and the hadronisation correction: compare [R ...]
compare() {
   cd "$here"
   for R in ${@:-0.2 0.3 0.4 0.5}; do
      root -l -b -q "compare.C+(\"$R\")"
   done
   root -l -b -q 'chad.C+'
}

case "$1" in
   build)
      build
      ;;
   bin)
      shift
      bin "$@"
      ;;
   gen)
      shift
      gen "$@"
      ;;
   chad)
      shift
      chad "$@"
      ;;
   spectra)
      shift
      spectra "$@"
      ;;
   compare)
      shift
      compare "$@"
      ;;
   all)
      for g in $GENS; do
         gen $g
      done
      chad
      spectra
      compare
      ;;
   *)
      sed -n 2,24p "$0"
      exit 1
      ;;
esac
