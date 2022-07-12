#ifndef RlBurninEstimationConvergenceAssessment_H
#define RlBurninEstimationConvergenceAssessment_H

#include "TypedDagNode.h"
#include "RevObject.h"
#include "WorkspaceObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev-wrapper of BurninEstimationConvergenceAssessment class.
     *
     * This class is the Rev-wrapper for BurninEstimationConvergenceAssessment core objetcs. 
     * See BurninEstimationConvergenceAssessment.h for more information.
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-09-06
     *
     */
    class BurninEstimationConvergenceAssessment : public WorkspaceObject {
        
    public:
        
        BurninEstimationConvergenceAssessment(void);                                                                                                //!< Default constructor
        
        // Basic utility functions
        virtual BurninEstimationConvergenceAssessment*      clone(void) const;                                                                      //!< Clone object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getConstructorFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                             getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                        printValue(std::ostream& o, bool user) const;                                                      //!< Print value (for user)
        
        // Member method inits
        RevPtr<RevVariable>                                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions
        
    protected:
        
        void                                                readTrace(const RevBayesCore::path &fn, std::vector<RevBayesCore::TraceNumeric> &data);
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable
        
        std::string                                         delimiter;
        std::set<RevBayesCore::path>                        filenames;
        std::string                                         burninMethod;
        bool                                                verbose;
    };
    
}

#endif
