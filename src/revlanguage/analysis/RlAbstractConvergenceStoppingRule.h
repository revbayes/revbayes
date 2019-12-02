#ifndef RlAbstractConvergenceStoppingRule_H
#define RlAbstractConvergenceStoppingRule_H

#include "BurninEstimatorContinuous.h"
#include "VariableMonitor.h"
#include "RlStoppingRule.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * RevLanguage abstract base class for convergence stopping rules.
     *
     *
     * This abstract class provides the base for (all) convergence stopping rules.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-09, version 1.0
     */
    class AbstractConvergenceStoppingRule : public StoppingRule {
        
    public:
        
        AbstractConvergenceStoppingRule(void);                                                                                           //!< Default constructor
        
        // Basic utility functions
        virtual StoppingRule*                       clone(void) const = 0;                                                              //!< Clone object
        virtual void                                constructInternalObject(void) = 0;                                                  //!< We construct the a new internal object.

        RevBayesCore::BurninEstimatorContinuous*    constructBurninEstimator(void);
        
        static const std::string&                   getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                 //!< Get class type spec
        virtual const MemberRules&                  getParameterRules(void) const;                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                //!< Get language type of the object
        
    protected:
        
        virtual void                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);     //!< Set member variable
        
        RevPtr<const RevVariable>                   burninMethod;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   frequency;
        
    };
    
}

#endif
