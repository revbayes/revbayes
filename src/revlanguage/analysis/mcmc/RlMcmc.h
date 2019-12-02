#ifndef RlMcmc_H
#define RlMcmc_H

#include "Mcmc.h"
#include "RlMonteCarloAnalysis.h"

#include <string>

namespace RevLanguage {
    
    
    /**
     * RevLanguage wrapper class for the power posterior analysis object.
     *
     *
     * The wraper class provides the Rev interface to the core class Mcmc.
     * See Mcmc.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class Mcmc : public MonteCarloAnalysis {
        
    public:
        
        Mcmc(void);                                                                                                               //!< Default constructor
        
        // Basic utility functions
        virtual Mcmc*                                   clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal PowerPosterior object.
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getConstructorFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        virtual const MemberRules&                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        
    protected:
        
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        virtual void                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable
        
        std::vector<std::string>                        getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                 //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                              //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                   //!< Get the title of this help entry
        
        RevPtr<const RevVariable>                       likelihood_heat;
        RevPtr<const RevVariable>                       posterior_heat;
        RevPtr<const RevVariable>                       prior_heat;
        
    };
    
}

#endif
