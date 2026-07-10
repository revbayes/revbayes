#ifndef RlPseudoDataLikelihood_H
#define RlPseudoDataLikelihood_H

#include <iostream>
#include <vector>

#include "PseudoDataLikelihood.h"
#include "ModelObject.h"
#include "TypedDagNode.h"


namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of PseudoDataLikelihood values.
     *
     * Currently (2026) the only purpose of this class is to store pseudodata likelihoods.
     * Therefore we do NOT want to allow interconverting between Real <-> PseudoDataLikelihood.
     *
     * @copyright Copyright 2026-
     * @author The RevBayes Development Core Team (Benjamin Redelings)
     * @since 2026-01-09
     *
     */
    class PseudoDataLikelihood : public ModelObject<RevBayesCore::PseudoDataLikelihood>  {
        
    public:
                                            PseudoDataLikelihood(void);                                                                       //!< Default constructor
                                            PseudoDataLikelihood(const RevBayesCore::PseudoDataLikelihood& l);
                                            PseudoDataLikelihood(RevBayesCore::PseudoDataLikelihood* l);
                                            PseudoDataLikelihood(RevBayesCore::TypedDagNode<RevBayesCore::PseudoDataLikelihood> *l);                    //!< Construct from DAG node
        
        // the value type definition
        PseudoDataLikelihood*                         clone(void) const override;                                                             //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        virtual const TypeSpec&             getTypeSpec(void) const;                                                                //!< Get language type of the object

        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
    };

}

#endif
