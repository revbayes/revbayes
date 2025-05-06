#ifndef RlPseudoData_H
#define RlPseudoData_H

#include "ModelObject.h"
#include "PseudoData.h"

#include <iostream>
#include <vector>

namespace RevLanguage {

    /**
     * @brief PseudoData: specifying data in terms of information it gives about a parameter via a likelihood function.
     */
    template <typename rlType>
    class PseudoData : public ModelObject<RevBayesCore::PseudoData<typename rlType::valueType> > {

    public:
        typedef typename RevBayesCore::PseudoData<typename rlType::valueType> valueType;

        PseudoData() {}
        PseudoData(const valueType& v): ModelObject<valueType>(v.clone()) {}
        PseudoData(valueType* v): ModelObject<valueType>(v) {}
        PseudoData(RevBayesCore::TypedDagNode<valueType>* n): ModelObject<valueType>(n) {}

        PseudoData<rlType>*         clone(void) const override
        {
            return new PseudoData<rlType>(*this);
        }

        static const std::string&   getClassType()
        {
            static std::string rev_type = "PseudoData<" + rlType::getClassType() + ">";

            return rev_type;
        }

        static const TypeSpec&      getClassTypeSpec(void)
        {
            static TypeSpec rev_type_spec = TypeSpec( getClassType(), &ModelObject<RevBayesCore::PseudoData<typename rlType::valueType> >::getClassTypeSpec(), &rlType::getClassTypeSpec() );

            return rev_type_spec;
        }

        const TypeSpec&             getTypeSpec(void) const override
        {
            return getClassTypeSpec();
        }


        // Member object functions
        std::string                 getGuiName(void) override { return "PseudoData"; }
        std::string                 getGuiUnicodeSymbol(void) override { return "PseudoData"; }
        std::string                 getGuiInfo(void) override { return ""; }
    };
}

#endif
