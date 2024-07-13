#ifndef Dist_InversePhylo_H
#define Dist_InversePhylo_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDistribution.h"
#include "Dist_phyloCTMC.h"

namespace RevLanguage {
    class Dist_InversePhylo : public TypedDistribution<AbstractHomologousDiscreteCharacterData> {
        
    public:
        Dist_InversePhylo(void);
        Dist_InversePhylo(const Dist_phyloCTMC& base_dist);
        virtual ~Dist_InversePhylo();
        
        // Basic utility functions
        Dist_InversePhylo*                               clone(void) const override;
        static const std::string&                        getClassType(void);
        static const TypeSpec&                           getClassTypeSpec(void);
        const TypeSpec&                                  getTypeSpec(void) const override;
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>*      createDistribution(void) const override;
        
        // Distribution-specific methods
        const MemberRules&                               getParameterRules(void) const override;
        std::vector<std::string>                         getDistributionFunctionAliases(void) const override;
        std::string                                      getDistributionFunctionName(void) const override;
        MethodTable                                      getDistributionMethods(void) const override;
        void                                             printValue(std::ostream& o) const override;
        
        // New method for calculating log probability
        double                                           calcLnProbability(void) const;
        
    protected:
        void                                             setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) override;
        
    private:
        Dist_phyloCTMC                                   base_distribution;
    };
    
}
#endif