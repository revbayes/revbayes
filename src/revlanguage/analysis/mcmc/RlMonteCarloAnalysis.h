#ifndef RlMonteCarloAnalysis_H
#define RlMonteCarloAnalysis_H

#include "MonteCarloAnalysis.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev-wrapper of MonteCarloAnalysis class.
     *
     * This class is the Rev-wrapper for MonteCarloAnalysis core objetcs. See MonteCarloAnalysis.h for more information.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-09-06
     *
     */
    class MonteCarloAnalysis : public WorkspaceToCoreWrapperObject<RevBayesCore::MonteCarloAnalysis> {
        
    public:
        
        MonteCarloAnalysis(void);                                                                                                               //!< Default constructor
        MonteCarloAnalysis(RevBayesCore::MonteCarloAnalysis *m);                                                                                //!< Constructor with internal core object
        
        // Basic utility functions
        virtual MonteCarloAnalysis*                 clone(void) const;                                                                          //!< Clone object
        static const std::string&                   getClassType(void);                                                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                     //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                              //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                    //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);     //!< Override to map member methods to internal functions
        
    protected:
        
        void                                        initializeMethods(void);                                                                    //!< Initialize the member methods
        virtual void                                printValue(std::ostream& o) const;                                                          //!< Print value
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);           //!< Set member variable
        
        RevPtr<const RevVariable>                   model;
        RevPtr<const RevVariable>                   moves;
        RevPtr<const RevVariable>                   monitors;
        RevPtr<const RevVariable>                   moveschedule;
        RevPtr<const RevVariable>                   num_init_attempts;
        RevPtr<const RevVariable>                   combine_traces;
        RevPtr<const RevVariable>                   num_runs;

        
    };
    
}

#endif
