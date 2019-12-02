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

        std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable
        
    private:
        RevPtr<const RevVariable>                       p;                                                                                      //!< The probabilities of the categories
        
    };
    
}

#endif
