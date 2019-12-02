#ifndef RlRateGeneratorSequence_H
#define RlRateGeneratorSequence_H

#include "ModelObject.h"
#include "RateGeneratorSequence.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the rate-map.
     *
     * The RevLanguage wrapper of the rate-map simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the RateGeneratorSequence.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-03-27, version 1.0
     *
     */
    class RateGeneratorSequence : public ModelObject<RevBayesCore::RateGeneratorSequence> {
        
    public:
        
        RateGeneratorSequence(void);                                                                                                              //!< Default constructor
        RateGeneratorSequence(RevBayesCore::RateGeneratorSequence *m);                                                                                          //!< Construct from RC::RateGeneratorSequence
        RateGeneratorSequence(RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence> *d);                                                              //!< Construct from DAG node
        
        // Basic utility functions
        RateGeneratorSequence*                            clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
    };
    
}

#endif
