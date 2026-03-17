#include "Dist_mkPrime.h"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Natural.h"
#include "OptionRule.h"
#include "PhyloCTMCSiteHomogeneousMkPrime.h"
#include "RbException.h"
#include "Real.h"
#include "RlString.h"
#include "RlTree.h"
#include "Tree.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Dist_mkPrime::Dist_mkPrime() : TypedDistribution<AbstractHomologousDiscreteCharacterData>() {
}

Dist_mkPrime::~Dist_mkPrime() {
}

Dist_mkPrime* Dist_mkPrime::clone(void) const {
    return new Dist_mkPrime(*this);
}

RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* Dist_mkPrime::createDistribution(void) const {
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree&>(tree->getRevObject()).getDagNode();
    size_t k_obs = static_cast<size_t>(static_cast<const Natural&>(kObs->getRevObject()).getValue());
    const RevBayesCore::TypedDagNode<std::int64_t>* u_node = static_cast<const Natural&>(u->getRevObject()).getDagNode();
    const std::string& prior_family = static_cast<const RlString&>(prior->getRevObject()).getValue();
    double prior_param = static_cast<const Real&>(priorParam->getRevObject()).getValue();
    size_t k_max = static_cast<size_t>(static_cast<const Natural&>(kMax->getRevObject()).getValue());
    const std::string& code = static_cast<const RlString&>(coding->getRevObject()).getValue();

    if (code != "all" && code != "variable" && code != "informative") {
        throw RbException("Invalid coding option \"" + code + "\" for dnMkPrime. Available codings: all, informative, variable.");
    }

    RevBayesCore::AscertainmentBias::Coding coding = RevBayesCore::AscertainmentBias::VARIABLE;
    if (code == "all") {
        coding = RevBayesCore::AscertainmentBias::ALL;
    } else if (code == "informative") {
        coding = RevBayesCore::AscertainmentBias::INFORMATIVE;
    }

    std::vector<double> prior_params(1, prior_param);

    return new RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime(
        tau,
        k_obs,
        u_node,
        prior_family,
        prior_params,
        k_max,
        coding);
}

const std::string& Dist_mkPrime::getClassType(void) {
    static std::string rev_type = "Dist_mkPrime";

    return rev_type;
}

const TypeSpec& Dist_mkPrime::getClassTypeSpec(void) {
    static TypeSpec rev_type_spec = TypeSpec(getClassType(), new TypeSpec(TypedDistribution<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec()));

    return rev_type_spec;
}

std::string Dist_mkPrime::getDistributionFunctionName(void) const {
    return "MkPrime";
}

const MemberRules& Dist_mkPrime::getParameterRules(void) const {
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if (!rules_set) {
        dist_member_rules.push_back(new ArgumentRule("tree", Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        dist_member_rules.push_back(new ArgumentRule("kObs", Natural::getClassTypeSpec(), "Observed number of states in this partition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY));
        dist_member_rules.push_back(new ArgumentRule("u", Natural::getClassTypeSpec(), "Number of unobserved states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
        dist_member_rules.push_back(new ArgumentRule("prior", RlString::getClassTypeSpec(), "Prior family on u.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("geometric")));
        dist_member_rules.push_back(new ArgumentRule("priorParam", Real::getClassTypeSpec(), "Prior parameter for the chosen prior family.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(0.4)));

        std::vector<std::string> coding_options;
        coding_options.push_back("all");
        coding_options.push_back("informative");
        coding_options.push_back("variable");
        dist_member_rules.push_back(new OptionRule("coding", new RlString("variable"), coding_options, "Ascertainment bias correction."));

        dist_member_rules.push_back(new ArgumentRule("kMax", Natural::getClassTypeSpec(), "Maximum allowed number of states.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(18)));

        rules_set = true;
    }

    return dist_member_rules;
}

const TypeSpec& Dist_mkPrime::getTypeSpec(void) const {
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}

void Dist_mkPrime::printValue(std::ostream& o) const {
    o << "MkPrime(tree=";
    o << (tree != NULL ? tree->getName() : "?");
    o << ", kObs=";
    o << (kObs != NULL ? kObs->getName() : "?");
    o << ", u=";
    o << (u != NULL ? u->getName() : "?");
    o << ", prior=";
    o << (prior != NULL ? prior->getName() : "?");
    o << ", priorParam=";
    o << (priorParam != NULL ? priorParam->getName() : "?");
    o << ", coding=";
    o << (coding != NULL ? coding->getName() : "?");
    o << ", kMax=";
    o << (kMax != NULL ? kMax->getName() : "?");
    o << ")";
}

void Dist_mkPrime::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var) {
    if (name == "tree") {
        tree = var;
    } else if (name == "kObs") {
        kObs = var;
    } else if (name == "u") {
        u = var;
    } else if (name == "prior") {
        prior = var;
    } else if (name == "priorParam") {
        priorParam = var;
    } else if (name == "coding") {
        coding = var;
    } else if (name == "kMax") {
        kMax = var;
    } else {
        Distribution::setConstParameter(name, var);
    }
}
