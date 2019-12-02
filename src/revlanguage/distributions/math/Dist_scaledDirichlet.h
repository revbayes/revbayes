#ifndef Dist_scaledDirichlet_H
#define Dist_scaledDirichlet_H

#include "ScaledDirichletDistribution.h"
#include "RlTypedDistribution.h"
#include "RlSimplex.h"

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the Scaled Dirichlet distribution.
     *
     * The RevLanguage wrapper of the Scaled Dirichlet distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the ScaledDirichletDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author Andrew Magee
     * @since 2017-03-23, version 1.0
     *
     */
    class Dist_scaledDirichlet : public TypedDistribution<Simplex> {
        
    public:
        Dist_scaledDirichlet( void );
        virtual ~Dist_scaledDirichlet();
        
        // Basic utility functions
        Dist_scaledDirichlet*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::ScaledDirichletDistribution*      createDistribution(void) const;
        
    protected:
        
        std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       beta;
        
    };
    
}

#endif
