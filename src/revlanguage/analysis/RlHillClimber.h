#ifndef RlHillClimber_H
#define RlHillClimber_H

#include "HillClimber.h"
#include "RlMaximumLikelihoodAnalysis.h"

#include <string>

namespace RevLanguage {
    
    
    /**
     * RevLanguage wrapper class for the maximum likelihood analysis object using a hill-climbing algorithm.
     *
     *
     * The wraper class provides the Rev interface to the core class HillClimber.
     * See HillClimber.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class HillClimber : public MaximumLikelihoodAnalysis {
        
    public:
        
        HillClimber(void);                                                                                                               //!< Default constructor
        
        // Basic utility functions
        virtual HillClimber*                            clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal PowerPosterior object.
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getConstructorFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        virtual const MemberRules&                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        
    protected:
        
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        virtual void                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable

    };
    
}

#endif
