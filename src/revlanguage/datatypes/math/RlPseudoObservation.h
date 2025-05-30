#ifndef RlPseudoObservation_H
#define RlPseudoObservation_H

#include "ModelObject.h"
#include "PseudoObservation.h"

#include <iostream>
#include <vector>

namespace RevLanguage {

    /**
     * @brief PseudoObservation: specifying data in terms of information it gives about a parameter via a likelihood function.
     */
    class PseudoObservation : public ModelObject<RevBayesCore::PseudoObservation> {

    public:
        typedef typename RevBayesCore::PseudoObservation valueType;

        PseudoObservation() {}
        PseudoObservation(const valueType& v): ModelObject<valueType>(v.clone()) {}
        PseudoObservation(valueType* v): ModelObject<valueType>(v) {}
        PseudoObservation(RevBayesCore::TypedDagNode<valueType>* n): ModelObject<valueType>(n) {}

        PseudoObservation*         clone(void) const override
        {
            return new PseudoObservation(*this);
        }

        static const std::string&   getClassType()
        {
            static std::string rev_type = "PseudoObservation";

            return rev_type;
        }

        static const TypeSpec&      getClassTypeSpec(void)
        {
            static TypeSpec rev_type_spec = TypeSpec( getClassType(), &ModelObject<RevBayesCore::PseudoObservation>::getClassTypeSpec() );

            return rev_type_spec;
        }

        const TypeSpec&             getTypeSpec(void) const override
        {
            return getClassTypeSpec();
        }


        // Member object functions
        std::string                 getGuiName(void) override { return "PseudoObservation"; }
        std::string                 getGuiUnicodeSymbol(void) override { return "PseudoObservation"; }
        std::string                 getGuiInfo(void) override { return ""; }
    };
}

#endif
