#ifndef Dist_BirthDeathJeffreys_H
#define Dist_BirthDeathJeffreys_H

#include "BirthDeathJeffreysDistribution.h"
#include "RlTypedDistribution.h"
#include "ModelVector.h"
#include "Real.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the multivariate normal distribution.
     *
     * The RevLanguage wrapper of the multivariate normal distribution takes care of create the internal distribution object
     * and provides the RevLanguage object that can be used within Rev.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck and Risa Takenaka)
     * @since 2012-08-06, version 1.0
     *
     */
    class Dist_BirthDeathJeffreys :  public TypedDistribution<ModelVector<RealPos> > {
    
    public:
                                                        Dist_BirthDeathJeffreys( void );
    
        // Basic utility functions
        Dist_BirthDeathJeffreys*                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
//        // Member method functions
//        virtual RevPtr<RevVariable>                     executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
//        virtual const MethodTable&                      getMethods(void) const;                                                                     //!< Get member methods

    
        // Distribution functions you have to override
        RevBayesCore::BirthDeathJeffreysDistribution*   createDistribution(void) const;
    
    protected:
    
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    
    private:

        RevPtr<const RevVariable>                       rho;                                                                                    //!< The taxon sampling fraction
        RevPtr<const RevVariable>                       condition;                                                                          //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                       start_age;                                                                           //!< The start time of the process since the origin
        std::string                                     start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                       limit;

    };
    
}

#endif
