#ifndef Dist_DirichletTimeTree_H
#define Dist_DirichletTimeTree_H

#include "DirichletTimeTreeDistribution.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev Dirichlet time tree distribution
     *
     * @brief Declaration of the Rev Dirichlet time tree distribution
     *
     * @author Fredrik Ronquist
     */
 class Dist_DirichletTimeTree : public TypedDistribution<TimeTree> {
        
    public:
                Dist_DirichletTimeTree( void );        //!< Constructor
        virtual ~Dist_DirichletTimeTree();             //!< Virtual destructor
        
        // Basic utility functions
        Dist_DirichletTimeTree*                         clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::DirichletTimeTreeDistribution*    createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       root_age;
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       taxa;
    };

}

#endif
