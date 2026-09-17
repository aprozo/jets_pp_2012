// variants.h — the systematic variations of the jet cross section, one table.
//
// Each variation is a named change of the inputs of the unfolding; the name is the tag of every
// derived file (data_blocks_R<R>_<name>.root, response_blocks_R<R>_e23_<name>.root, xsec_..._<name>).
// The components and their members follow the published analysis:
//   energy scale   response only: detector-jet pT shifted by +-pT sqrt(((1-R_T) 0.011)^2 + (R_T 0.032)^2)
//                  and/or the jet-patch (and high-tower) DSM thresholds moved by +-1 ADC in the jet-to-patch
//                  match (the level keeps the nominal simulator decision); six members, envelope
//   tracking       data only: 1 % of the tracks of every jet removed, i.e. pT x (1 - 0.01 (1 - R_T));
//                  one member, one-sided envelope
//   underlying     detector-level UE subtraction x 0.86 / 1.18 in data and response, particle level
//   event          untouched; the seven (data, response) combinations of the published analysis, envelope
//   embedding      1000 Poisson replicas of the entry counts of the fine response, fakes and misses,
//   statistics     each times the cell's average weight, re-unfolded; per-bin RMS, symmetric
//   regularisation Bayesian only: iterations 2 and 6 around the nominal 4, and the prior tilted by
//                  (pT / 20 GeV)^(+-0.5); one envelope
// Up and down deviations from the nominal are combined in quadrature separately.
#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace Syst {

// the size of each variation, in one place
const double kJesNeutralTerm = 0.011; // relative tower-energy scale uncertainty
const double kJesChargedTerm = 0.032; // relative track-momentum scale uncertainty
const double kTrackLoss = 0.01;       // fraction of the tracks of a jet removed
const double kUeDown = 0.86;          // UE subtraction scaled down
const double kUeUp = 1.18;            // UE subtraction scaled up
const double kPriorTilt = 0.5;        // exponent of the (pT / 20 GeV) prior tilt

struct Variant {
   std::string name;       // "" = nominal
   double jesSign = 0;     // +1 / -1: detector-jet pT shift (response)
   int thrShift = 0;       // DSM threshold shift in ADC counts (response)
   double ueResp = 1.0;    // detector-level UE subtraction factor, response
   double ueData = 1.0;    // detector-level UE subtraction factor, data
   double trkThin = 0.0;   // fraction of tracks removed from the data jets
   double priorTilt = 0.0; // truth (prior) reweighted by (pT/20)^priorTilt
   bool embStat = false;   // replica loop of the embedding statistics
   bool ueParticle = true; // particle-level UE subtraction (off only for the no-UE jet definition)
   std::string rtag;       // tag of the response file this variant needs (set by Parse)
   std::string dtag;       // tag of the data file this variant needs (set by Parse)
   std::string respTag() const { return rtag; }
   std::string dataTag() const { return dtag; }
};

// the detector-jet pT shift of the energy-scale members, as a fraction of the raw jet pT
inline double JesShift(double neutralFraction, double sign)
{
   return sign * std::sqrt(std::pow((1 - neutralFraction) * kJesNeutralTerm, 2) +
                           std::pow(neutralFraction * kJesChargedTerm, 2));
}

// which response and data files a variant reads: a variation that changes the response gets a response
// tag, one that changes the data a data tag, and the nominal neither
inline void SetFileTags(Variant &variant)
{
   if (variant.jesSign > 0) variant.rtag += "_jesUp";
   else if (variant.jesSign < 0) variant.rtag += "_jesDn";
   if (variant.thrShift > 0) variant.rtag += "_thrP1";
   else if (variant.thrShift < 0) variant.rtag += "_thrM1";
   if (variant.ueResp != 1.0) variant.rtag += variant.ueResp < 1.0 ? "_ueR086" : "_ueR118";
   if (variant.trkThin != 0.0) variant.dtag += "_trk1";
   if (variant.ueData != 1.0) variant.dtag += variant.ueData < 1.0 ? "_ueD086" : "_ueD118";
}

// one member of the band, by name (the nominal is "")
inline Variant ParseMember(const std::string &name)
{
   Variant variant;
   variant.name = name;
   if (name.empty()) {
      return variant;
   } else if (name == "jesUp") {
      variant.jesSign = +1;
   } else if (name == "jesDn") {
      variant.jesSign = -1;
   } else if (name == "thrP1") {
      variant.thrShift = +1;
   } else if (name == "thrM1") {
      variant.thrShift = -1;
   } else if (name == "jesUp_thrP1") {
      variant.jesSign = +1;
      variant.thrShift = +1;
   } else if (name == "jesDn_thrM1") {
      variant.jesSign = -1;
      variant.thrShift = -1;
   } else if (name == "trk1") {
      variant.trkThin = kTrackLoss;
   } else if (name == "ueD086") {
      variant.ueData = kUeDown;
   } else if (name == "ueD118") {
      variant.ueData = kUeUp;
   } else if (name == "ueR086") {
      variant.ueResp = kUeDown;
   } else if (name == "ueR118") {
      variant.ueResp = kUeUp;
   } else if (name == "ueD086_ueR086") {
      variant.ueData = kUeDown;
      variant.ueResp = kUeDown;
   } else if (name == "ueD118_ueR118") {
      variant.ueData = kUeUp;
      variant.ueResp = kUeUp;
   } else if (name == "tiltUp") {
      variant.priorTilt = +kPriorTilt;
   } else if (name == "tiltDn") {
      variant.priorTilt = -kPriorTilt;
   } else if (name == "embstat") {
      variant.embStat = true;
   } else {
      throw std::runtime_error("unknown variant " + name);
   }
   SetFileTags(variant);
   return variant;
}

// Every variant by name: a member, the no-UE jet definition "noue", or "noue_<member>" (that member on
// that definition). On the no-UE definition the UE members add or remove the same 14 / 18 % of the
// underlying event to the raw jets instead of scaling a subtraction that is not there.
inline Variant Parse(const std::string &name)
{
   const bool noUe = name == "noue" || name.rfind("noue_", 0) == 0;
   const std::string memberName = noUe ? (name == "noue" ? "" : name.substr(5)) : name;
   Variant member = ParseMember(memberName);
   if (!noUe) return member;
   Variant variant = member;
   variant.name = name;
   variant.ueParticle = false;
   variant.ueData = member.ueData == 1.0 ? 0.0 : member.ueData - 1.0;
   variant.ueResp = member.ueResp == 1.0 ? 0.0 : member.ueResp - 1.0;
   variant.rtag = "_noue" + member.rtag;
   variant.dtag = "_noue" + member.dtag;
   return variant;
}

// the response variants to build (one build_resp pass each) and the data variants (one pass for all)
inline std::vector<std::string> ResponseVariants()
{
   return {"jesUp", "jesDn", "thrP1", "thrM1", "jesUp_thrP1", "jesDn_thrM1", "ueR086", "ueR118"};
}

inline std::vector<std::string> DataVariants()
{
   return {"trk1", "ueD086", "ueD118"};
}

// alternative jet definitions (not members of the band); the data pass also builds every data variant on them
inline std::vector<std::string> DefinitionVariants()
{
   return {"noue"};
}

// every data build of one pass: the nominal, the data variants, and both on each jet definition
inline std::vector<std::string> DataBuilds()
{
   std::vector<std::string> builds = {""};
   for (const std::string &data : DataVariants()) builds.push_back(data);
   for (const std::string &definition : DefinitionVariants()) {
      builds.push_back(definition);
      for (const std::string &data : DataVariants()) builds.push_back(definition + "_" + data);
   }
   return builds;
}

// one component of the band: its name, its members (variant names; "" = nominal) and whether the members
// give an envelope or a symmetric RMS. The Bayes iterations are handled by the driver (nominal 4, the
// members 2 and 6 in the regularisation component).
struct Component {
   std::string name;                 // the label of the component in the table
   std::vector<std::string> members; // the variant names it is built from
   bool rms;                         // true: the member carries its deviation as a per-bin RMS

   Component(const std::string &componentName, const std::vector<std::string> &memberNames, bool isRms = false)
      : name(componentName), members(memberNames), rms(isRms)
   {
   }
};

inline std::vector<Component> Components(bool bayes)
{
   std::vector<Component> components = {
      {"energy scale", {"jesDn_thrM1", "jesUp_thrP1", "jesDn", "jesUp", "thrM1", "thrP1"}},
      {"tracking", {"trk1"}},
      {"embedding statistics", {"embstat"}, true},
      {"underlying event", {"ueD086", "ueR086", "ueD086_ueR086", "ueR118", "ueD118", "ueD118_ueR118"}},
   };
   if (bayes) components.push_back({"regularisation", {"iter2", "iter6", "tiltUp", "tiltDn"}});
   return components;
}

} // namespace Syst
