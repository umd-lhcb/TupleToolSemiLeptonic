// Include files

// from Gaudi
#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"
#include "GaudiKernel/SmartIF.h"

// from LHCb
#include "Event/Particle.h"
#include "Event/RecVertex.h"
#include "Event/Vertex.h"

// from Phys
#include "Kernel/GetIDVAlgorithm.h"
#include "Kernel/IBTaggingTool.h"
#include "Kernel/IDVAlgorithm.h"

// local
#include "TupleToolTagging.h"

using namespace LHCb;

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_TOOL_FACTORY(TupleToolTagging)

//==============================================================================
// Standard constructor
//==============================================================================
TupleToolTagging::TupleToolTagging(const std::string& type,
                                   const std::string& name,
                                   const IInterface* parent)
    : TupleToolBase(type, name, parent) {
  declareInterface<IParticleTupleTool>(this);
}

//==============================================================================
// Initialize method
//==============================================================================
StatusCode TupleToolTagging::initialize() {
  if (!TupleToolBase::initialize()) return StatusCode::FAILURE;

  // get parent DaVinci algorithm
  m_parentDVA = Gaudi::Utils::getIDVAlgorithm(contextSvc(), this);
  if (m_parentDVA == nullptr) {
    return Error("Unable to retrieve parent DVAlg", StatusCode::FAILURE);
  }

  // write out information from FlavourTag saved on DST
  if (m_useFTfromDST) {
    info() << "Read tagging results saved in TES. " << endmsg;

    // sanity checks and warning for user
    if (m_addMVAFeatureInfo) {
      warning() << "Asked to read tagging results from TES. Will ignore "
                   "AddMVAFeatureInfo."
                << endmsg;
    }
    if (m_addTagPartsInfo) {
      warning() << "Asked to read tagging results from TES. Will ignore "
                   "AddTagPartsInfo."
                << endmsg;
    }
  }
  // configured to run FT, thus have to get/create TaggingTool
  else {
    // if no name is given, get TaggingTool from parent DVA, else create a new
    // private tool
    if (m_nameTaggingTool == "") {
      m_taggingTool = m_parentDVA->flavourTagging();
    } else {
      m_taggingTool = tool<IBTaggingTool>(m_nameTaggingTool, this);
    }
    if (m_taggingTool == nullptr) {
      return Error("Unable to retrieve the IBTaggingTool tool",
                   StatusCode::FAILURE);
    }

    // set active tagger list to the list of active taggers of the TaggingTool
    info() << "Will use TaggingTool " << m_taggingTool->name()
           << " for tagging." << endmsg;
    m_activeTaggerTypeNames = m_taggingTool->activeTaggerTypeNames();
  }

  // transform TaggerType names to pair of TaggerTypes and names
  info() << "Will save information for following tags (if available): ";
  std::vector<std::string> invalidTaggerNames;
  for (auto activeTaggerName : m_activeTaggerTypeNames) {
    auto activeTaggerType = Tagger::TaggerTypeToType(activeTaggerName);
    if (activeTaggerType != Tagger::TaggerType::unknown &&
        activeTaggerType != Tagger::TaggerType::none) {
      info() << activeTaggerName << " ";
      m_activeTaggers.emplace(activeTaggerName, activeTaggerType);
    } else {
      invalidTaggerNames.push_back(activeTaggerName);
    }
  }
  info() << "." << endmsg;

  if (!(invalidTaggerNames.empty())) {
    warning() << "." << endmsg;
    warning() << "Could not add the following taggers, as they are unknown:";
    for (auto invalidTaggerName : invalidTaggerNames) {
      warning() << " " << invalidTaggerName;
    }
    warning() << "." << endmsg;
  }

  if (m_tagNonBeauty) {
    info() << "Running Flavour Tagging on all particles (not only b-hadrons)"
           << endmsg;
  }

  return StatusCode::SUCCESS;
}

//==============================================================================
// fill method
//==============================================================================
StatusCode TupleToolTagging::fill(const Particle* motherPart,
                                  const Particle* signalPart,
                                  const std::string& head,
                                  Tuples::Tuple& tuple) {
  // consistency checks
  Assert(signalPart && motherPart && m_parentDVA && m_taggingTool,
         "Should not happen, you are inside TupleToolTagging.cpp");

  // nothing to tag on something which is not a B (unless explicitly required)
  if (!m_tagNonBeauty) {
    if (!(signalPart->particleID().hasBottom())) return StatusCode::SUCCESS;
  }

  // instantiate empty FlavourTag and FlavourTags
  FlavourTag flavourTag;
  FlavourTags* flavourTags = nullptr;

  // create successful statuscode
  StatusCode sc = StatusCode::SUCCESS;
  bool check = false;

  // if requested, retrieve FlavourTag objects from DST location
  if (m_useFTfromDST) {
    info() << "Will read flavour tagging results from TES." << endmsg;
    // get location of signal particles, then replace with expected location of
    // FlavourTags vector
    std::string signalPartLoc = objectLocation(signalPart->parent());
    boost::replace_all(signalPartLoc, "/Particles", "/FlavourTags");

    // get FlavourTags from expected DST location
    if (exist<LHCb::FlavourTags>(signalPartLoc, IgnoreRootInTES)) {
      flavourTags = get<LHCb::FlavourTags>(signalPartLoc, IgnoreRootInTES);
    }
  }

  // check if retrieved FlavourTags belong to to the signal particle
  if (flavourTags != nullptr) {
    for (FlavourTags::const_iterator it = flavourTags->begin();
         it != flavourTags->end(); ++it) {
      if (signalPart != (**it).taggedB()) continue;
      flavourTag = **it;
      check = true;
    }
    if (!check) sc = StatusCode::FAILURE;
  } else {
    const RecVertex* assocVtx =
        dynamic_cast<const RecVertex*>(m_parentDVA->bestVertex(motherPart));
    if (assocVtx == nullptr) {
      sc = m_taggingTool->tag(flavourTag, signalPart);
    } else {
      sc = m_taggingTool->tag(flavourTag, signalPart, assocVtx);
    }
  }

  // replace if (!(sc.isSuccess())) return sc; to avoid situation where tool
  // returns FALIURE when bestPV is not found
  if (!(sc.isSuccess())) Warning("The tagging algorithm failed");
  // Start filling the tuples
  // get tuple branch prefix
  const std::string prefix = fullName(head);
  sc = fillTagInfo(flavourTag, prefix, tuple);

  return sc;
}

StatusCode TupleToolTagging::fillTagInfo(const FlavourTag& flavourTag,
                                         const std::string& prefix,
                                         Tuples::Tuple& tuple) {
  // retrieve taggers
  auto taggers = flavourTag.taggers();

  bool test = true;
  std::string tuplePrefix = prefix;
  std::string taggerName("none");
  Tagger::TaggerType taggerType(Tagger::TaggerType::none);
  // loop over map of all active taggers
  for (auto activeTaggerNameAndType : m_activeTaggers) {
    taggerName = activeTaggerNameAndType.first;
    taggerType = activeTaggerNameAndType.second;
    tuplePrefix = prefix + "_" + taggerName;

    // check if active tagger has a Tagger object in Taggers
    // loop over all Tagger object and check if type matches active tagger type
    auto tagger = std::find_if(taggers.begin(), taggers.end(),
                               [&taggerType](Tagger current_tagger) -> bool {
                                 return (current_tagger.type() == taggerType);
                               });
    bool hasTagged = (tagger != taggers.end());

    int taggerDec = 0;
    double taggerEta = 0.5;
    double taggerCharge = 0.;
    double taggerMVAout = -1.;
    if (hasTagged) {
      taggerDec = static_cast<int>(tagger->decision());
      taggerEta = static_cast<double>(tagger->omega());
    }

    // tagging decision and estimated mistag probability
    test &= tuple->column(tuplePrefix + "_TAGDEC", taggerDec);
    test &= tuple->column(tuplePrefix + "_TAGETA", taggerEta);

    // charge used for decision and raw (uncalibrated) MVA classifier output
    if (isVerbose()) {
      if (hasTagged) {
        taggerCharge = static_cast<double>(tagger->charge());
        taggerMVAout = static_cast<double>(tagger->mvaValue());
      }
      test &= tuple->column(tuplePrefix + "_CHARGE", taggerCharge);
      test &= tuple->column(tuplePrefix + "_MVAOUT", taggerMVAout);
    }

    if (!m_useFTfromDST) {
      if (m_addMVAFeatureInfo) {
        auto featureNames = m_taggingTool->featureNames(taggerType);
        std::vector<double> featureValues;
        if (hasTagged) {
          featureValues = std::move(m_taggingTool->featureValues(taggerType));
        } else {
          featureValues =
              std::move(std::vector<double>(featureNames.size(), -99999.));
        }
        // write features to tuple
        for (size_t i = 0; i < featureNames.size(); ++i) {
          test &=
              tuple->column(tuplePrefix + "_MVAFeature_" + featureNames.at(i),
                            featureValues.at(i));
        }
      }

      if (m_addTagPartsInfo) {
        auto featureNames = m_taggingTool->featureNamesTagParts(taggerType);
        std::vector<std::vector<double>> featureValueMatrix;
        if (hasTagged) {
          featureValueMatrix =
              std::move(m_taggingTool->featureValuesTagParts(taggerType));
        } else {
          featureValueMatrix = std::move(
              std::vector<std::vector<double>>(featureNames.size(), {-99999.}));
        }

        // write features to tuple
        const std::string tupleNumName =
            tuplePrefix + "_TagPartsFeature_" + "NUM";
        for (size_t i = 0; i < featureNames.size(); ++i) {
          test &= tuple->farray(
              tuplePrefix + "_TagPartsFeature_" + featureNames.at(i),
              featureValueMatrix.at(i), tupleNumName, m_maxTagPartsInfo);
        }
      }
    }
  }

  return StatusCode(test);
}
