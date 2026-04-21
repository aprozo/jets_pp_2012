TRIGGERS=(JP2 HT2)
# RADII=(0.2 0.3 0.4 0.5 0.6)
RADII=(0.5 0.6)

for TRIGGER in "${TRIGGERS[@]}"; do
  for R in "${RADII[@]}"; do
    echo "Unfolding for trigger: $TRIGGER and R=$R"
    filename=/home/prozorov/dev/star/jets_pp_2012/output/merged_matching_"$TRIGGER"_R"$R.root"
    root -l -b -q -e 'gSystem->Load("libRooUnfold");' "unfold.cxx+(\"$filename\")"
  done
done