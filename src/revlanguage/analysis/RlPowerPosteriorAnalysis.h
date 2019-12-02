#ifndef RlPowerPosteriorAnalysis_H
#define RlPowerPosteriorAnalysis_H

#include "PowerPosteriorAnalysis.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {

    
    /**
     * RevLanguage wrapper class for the power posterior analysis object.
     *
     *
     * The wraper class provides the Rev interface to the core class PowerPosteriorAnalysis.
     * See PowerPosteriorAnalysis.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class PowerPosteriorAnalysis : public WorkspaceToCoreWrapperObject<RevBayesCore::PowerPosteriorAnalysis> {

    public:

        PowerPosteriorAnalysis(void);                                                                                                               //!< Default constructor

        // Basic utility functions
        virtual PowerPosteriorAnalysis*             clone(void) const;                                                                      //!< Clone object
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
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable

        RevPtr<const RevVariable>                   model;
        RevPtr<const RevVariable>                   moves;
        RevPtr<const RevVariable>                   monitors;
        RevPtr<const RevVariable>                   powers;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   cats;
        RevPtr<const RevVariable>                   alphaVal;
        RevPtr<const RevVariable>                   sampFreq;
        RevPtr<const RevVariable>                   proc_per_lik;

    };

}

#endif
