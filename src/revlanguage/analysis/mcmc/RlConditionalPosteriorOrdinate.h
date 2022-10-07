#ifndef RlConditionalPosteriorOrdinate_H
#define RlConditionalPosteriorOrdinate_H

#include "ConditionalPosteriorOrdinate.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * @brief Rev wrapper of the Stepping-Stone-Sampler class.
     *
     * The stepping stone sampler analyzes the output of a power posterior and
     * computes the marginal likelihood.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class ConditionalPosteriorOrdinate : public WorkspaceToCoreWrapperObject<RevBayesCore::ConditionalPosteriorOrdinate> {
        
    public:
        
        ConditionalPosteriorOrdinate(void);                                                                                                 //!< Default constructor
        
        // Basic utility functions
        virtual ConditionalPosteriorOrdinate*       clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal PowerPosterior object.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getConstructorFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions
        
    protected:
        
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   delimiter;
        RevPtr<const RevVariable>                   col_names_to_skip;
        RevPtr<const RevVariable>                   filename;
        
    };
    
}

#endif
