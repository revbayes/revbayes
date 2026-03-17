#ifndef Dist_mkPrime_H
#define Dist_mkPrime_H

#include <iosfwd>
#include <string>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
class TypeSpec;

class Dist_mkPrime : public TypedDistribution<AbstractHomologousDiscreteCharacterData> {
public:
    Dist_mkPrime(void);
    virtual ~Dist_mkPrime();

    Dist_mkPrime* clone(void) const;
    static const std::string& getClassType(void);
    static const TypeSpec& getClassTypeSpec(void);
    std::string getDistributionFunctionName(void) const;
    const TypeSpec& getTypeSpec(void) const;
    const MemberRules& getParameterRules(void) const;
    void printValue(std::ostream& o) const;

    RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* createDistribution(void) const;

protected:
    void setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);

private:
    RevPtr<const RevVariable> tree;
    RevPtr<const RevVariable> kObs;
    RevPtr<const RevVariable> u;
    RevPtr<const RevVariable> prior;
    RevPtr<const RevVariable> priorParam;
    RevPtr<const RevVariable> coding;
    RevPtr<const RevVariable> kMax;
};

}

#endif
