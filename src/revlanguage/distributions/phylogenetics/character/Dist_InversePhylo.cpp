#include "Dist_InversePhylo.h"
#include "Dist_phyloCTMC.h"
#include "InversePhyloCTMC.h"
#include "HomologousDiscreteCharacterData.h" // Include the concrete class

using namespace RevLanguage;


const TypeSpec& Dist_InversePhylo::getTypeSpec(void) const {
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


template<class charType>
RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* 
Dist_InversePhylo<charType>::createDistribution(void) const {
    RevBayesCore::TypedDistribution<RevBayesCore::HomologousDiscreteCharacterData<charType>>* base_dist = 
        dynamic_cast<RevBayesCore::TypedDistribution<RevBayesCore::HomologousDiscreteCharacterData<charType>>*>(base_distribution.createDistribution());
    if (base_dist == nullptr) {
        throw Exception("Failed to cast base distribution to correct type");
    }
    return new RevBayesCore::InversePhyloCTMC<RevBayesCore::HomologousDiscreteCharacterData<charType>>(*base_dist);
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
