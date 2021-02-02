#ifndef Dist_categorical_H
#define Dist_categorical_H

#include "CategoricalDistribution.h"
#include "Natural.h"
#include "Probability.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * \brief Categorical distribution
     *
     * This class is the language wrapper for the categorical distribution.
     * The distribution takes a simplex parameter, which defines the
     * probabilities of the categories, and the sample space. If the
     * simplex has k elements, then the sample space is { 1, 2, ..., k }.
     *
     * See RevBayesCore::CategoricalDistribution for more information.
     *
     * \author Fredrik Ronquist
     *
     */
    class Dist_categorical :  public TypedDistribution<Natural> {
        
    public:
        Dist_categorical( void );                                                                                                               //!< Default constructor
        virtual                                         ~Dist_categorical();                                                                    //!< Virtual destructor
        
        // Basic utility functions
        Dist_categorical*                               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type (static)
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get Rev type spec (static)
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get Rev type spec (dynamic)
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')

        // Distribution function overridden here
        RevBayesCore::CategoricalDistribution*          createDistribution(void) const;
        
    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable
        
    private:
        RevPtr<const RevVariable>                       p;                                                                                      //!< The probabilities of the categories
        
    };
    
}

#endif
