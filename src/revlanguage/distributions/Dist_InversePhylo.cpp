#include "Dist_InversePhylo.h"
#include "InversePhyloCTMCDistribution.h"

using namespace RevLanguage;

Dist_InversePhylo::Dist_InversePhylo(void) : TypedDistribution<AbstractHomologousDiscreteCharacterData>() {
}

Dist_InversePhylo::Dist_InversePhylo(const Dist_phyloCTMC& base_dist) : TypedDistribution<AbstractHomologousDiscreteCharacterData>(), base_distribution(base_dist) {
}

Dist_InversePhylo::~Dist_InversePhylo() {
}

Dist_InversePhylo* Dist_InversePhylo::clone(void) const {
    return new Dist_InversePhylo(*this);
}

const std::string& Dist_InversePhylo::getClassType(void) {
    static std::string rev_type = "Dist_InversePhylo";
    return rev_type;
}

const TypeSpec& Dist_InversePhylo::getClassTypeSpec(void) {
    static TypeSpec rev_type_spec = TypeSpec(getClassType(), new TypeSpec(TypedDistribution<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec()));
    return rev_type_spec;
}

const TypeSpec& Dist_InversePhylo::getTypeSpec(void) const {
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* Dist_InversePhylo::createDistribution(void) const {
    RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* base_dist = base_distribution.createDistribution();
    return new RevBayesCore::InversePhyloCTMCDistribution(*base_dist);
}

const MemberRules& Dist_InversePhylo::getParameterRules(void) const {
    return base_distribution.getParameterRules();
}

std::vector<std::string> Dist_InversePhylo::getDistributionFunctionAliases(void) const {
    std::vector<std::string> aliases = base_distribution.getDistributionFunctionAliases();
    aliases.push_back("invPhylo");
    return aliases;
}

std::string Dist_InversePhylo::getDistributionFunctionName(void) const {
    return "inversePhylo";
}

MethodTable Dist_InversePhylo::getDistributionMethods(void) const {
    return base_distribution.getDistributionMethods();
}

void Dist_InversePhylo::printValue(std::ostream& o) const {
    o << "InversePhylo(";
    base_distribution.printValue(o);
    o << ")";
}

void Dist_InversePhylo::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    base_distribution.setConstParameter(name, var);
}

double Dist_InversePhylo::calcLnProbability(void) const {
    return -base_distribution.calcLnProbability();
}
