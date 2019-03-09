#ifndef ANALYSIS_PHYS_DECAYTREETUPLE_TUPLETOOLTAGGING_H
#define ANALYSIS_PHYS_DECAYTREETUPLE_TUPLETOOLTAGGING_H 1

// Include files
// from Phys
#include "Kernel/IParticleTupleTool.h"  // Interface
#include "Kernel/IVertexFit.h"

// from Analysis
#include "DecayTreeTupleBase/TupleToolBase.h"

class IDVAlgorithm;
class IDistanceCalculator;
class IBTaggingTool;
class IVertexFit;

struct VerboseData {
  double id, p, px, py, pz, pt, theta, phi;
  double pid_e, pid_mu, pid_k, pid_p;
  double ip, chi2, bip, bchi2, bp_chi2;
  VerboseData()
      : id(0),
        p(0),
        px(0),
        pz(0),
        pt(0),
        theta(0),
        phi(0),
        pid_e(0),
        pid_mu(0),
        pid_k(0),
        pid_p(0),
        ip(0),
        chi2(0),
        bip(0),
        bchi2(0),
        bp_chi2(0) {}
};

/** @class TupleToolTagging TupleToolTagging.h
 *
 * \brief Fill the tagging information in the tuple.
 *
 * The output of this DecayTreeTupleTool consists of the tag decision (TAGDEC)
 * and the estimated mistag probability (TAGETA) of each tagging algorithm. The
 * list of considered tagging algorithms depends on whether the tag results are
 * read out from the TES or if they are produced by the associated BTaggingTool.
 *
 * If the results are read out from TES (UseFTfromDST is set to true), the list
 * of results written to the tuple depends on the list of active taggers defined
 * in the property ActiveTaggers. Elsewise, this list will be overwritten by the
 * list of active taggers in the associated BTaggingTool.
 *
 * If Verbose is set the effective charge of the tagging object and the MVA
 * output are filled for each tagger.
 *
 * @see  BTaggingTool
 * @see  Tagger
 *
 * @author Jeremie Borel
 * @author Julian Wishahi
 * @date   2007-2017
 */

class TupleToolTagging : public TupleToolBase,
                         virtual public IParticleTupleTool {
 public:
  /// Standard constructor
  TupleToolTagging(const std::string& type, const std::string& name,
                   const IInterface* parent);

  StatusCode initialize() override;

  StatusCode fill(const LHCb::Particle* motherPart,
                  const LHCb::Particle* signalPart, const std::string& head,
                  Tuples::Tuple& tuple) override;

 private:
  StatusCode fillTagInfo(const LHCb::FlavourTag& flavourTag,
                         const std::string& prefix, Tuples::Tuple& tuple);

  Gaudi::Property<std::string> m_nameTaggingTool{
      this, "TaggingToolName", "BTaggingTool/BTaggingTool",
      "Name of the TaggingTool to be used"};

  Gaudi::Property<bool> m_useFTfromDST{
      this, "UseFTfromDST", false,
      "Do not rerun flavour tagging but retrieve from DST"};

  Gaudi::Property<std::vector<std::string>> m_activeTaggerTypeNames{
      this,
      "ActiveTaggerTypes",
      {"OS_Muon", "OS_Electron", "OS_Kaon", "VtxCharge", "OS_Charm",
       "SS_nnetKaon", "SS_PionBDT", "SS_Proton"},
      "List of active tagger types. This list is only used when UseFTfromDST "
      "is set."};

  std::map<std::string, LHCb::Tagger::TaggerType> m_activeTaggers;

  Gaudi::Property<bool> m_tagNonBeauty{
      this, "TagNonBeauty", false,
      "Also run tagging on candidates that do not contain a beauty quark"};

  Gaudi::Property<bool> m_addMVAFeatureInfo{
      this, "AddMVAFeatureInfo", false,
      "Add feature values used by the MVA for the tagging decision to the "
      "tuple"};

  Gaudi::Property<bool> m_addTagPartsInfo{
      this, "AddTagPartsInfo", false,
      "Add feature values of all tagging particles considered"};

  Gaudi::Property<int> m_maxTagPartsInfo{
      this, "MaximumTagPartsInfo", 1000,
      "Maximum of rows to be saved when adding VerboseFeatureInfo"};

  IDVAlgorithm* m_parentDVA = nullptr;  ///< pointer to the parent DV algorithm
  IBTaggingTool* m_taggingTool = nullptr;  ///< pointer to the BTaggingTool
};

#endif  // ANALYSIS_PHYS_DECAYTREETUPLE_TUPLETOOLTAGGING_H
