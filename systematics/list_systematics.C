// list_systematics.C — print the committed systematic matrix
// (config.h::Systematics()) for the shell driver: one "SYST <name>
// <needsResponse>" line per variation. The single bridge from the config-as-code
// presets to run_systematics.sh, so the list lives in exactly one place.
// Usage: root -l -b -q list_systematics.C
#include "../new_ana/config.h"
using namespace CrossSectionConfig;

void list_systematics()
{
   for (const auto &s : Systematics())
      printf("SYST %s %d\n", s.name.c_str(), s.needsResponse() ? 1 : 0);
}
