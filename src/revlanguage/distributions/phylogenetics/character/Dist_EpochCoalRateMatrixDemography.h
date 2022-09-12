#ifndef Dist_EpochCoalRateMatrixDemography_H
#define Dist_EpochCoalRateMatrixDemography_H

#include "EpochCoalRateMatrixDemography.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the EpochCoalRateMatrixDemography distribution.
     *
     * The RevLanguage wrapper of the EpochCoalRateMatrixDemography distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the EpochCoalRateMatrixDemographyDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2022-09-12, version 1.2
     *
     */
    class Dist_EpochCoalRateMatrixDemography : public TypedDistribution< ModelVector<RealPos> > {
        
    public:
        Dist_EpochCoalRateMatrixDemography(void);
        virtual                                        ~Dist_EpochCoalRateMatrixDemography(void);
        
        // Basic utility functions
        Dist_EpochCoalRateMatrixDemography*             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::EpochCoalRateMatrixDemography*    createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       Ne;
        RevPtr<const RevVariable>                       mu;
        RevPtr<const RevVariable>                       times;
        RevPtr<const RevVariable>                       ancestral_sfs;
        RevPtr<const RevVariable>                       num_sites;
        RevPtr<const RevVariable>                       num_individuals;
        RevPtr<const RevVariable>                       folded;
        RevPtr<const RevVariable>                       coding;
        
    };
    
}

#endif
