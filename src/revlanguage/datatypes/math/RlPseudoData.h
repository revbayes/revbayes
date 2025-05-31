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

        PseudoData() { initMethods(); }
        PseudoData(const valueType& v): ModelObject<valueType>(v.clone()) { initMethods(); }
        PseudoData(valueType* v): ModelObject<valueType>(v) { initMethods(); }
        PseudoData(RevBayesCore::TypedDagNode<valueType>* n): ModelObject<valueType>(n) { initMethods(); }

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


        RevLanguage::RevPtr<RevLanguage::RevVariable> executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
        {
            if (name == "lnLikelihood")
            {
                found = true;
        
                double x = static_cast<const Real&>( args[0].getVariable()->getRevObject() ).getValue();

                double l = this->dag_node->getValue()(x);
                return RevPtr<RevVariable>( new RevVariable( new Real( l ) ) );
            }

            return ModelObject<RevBayesCore::PseudoData<typename rlType::valueType> >::executeMethod( name, args, found );
        }

        void initMethods( void )
        {
            ArgumentRules* lnLikelihood_rules = new ArgumentRules();
            lnLikelihood_rules->push_back( new ArgumentRule( "x", rlType::getClassTypeSpec(), "The likelihood for a parameter value.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
            this->methods.addFunction( new MemberProcedure( "lnLikelihood", Real::getClassTypeSpec(), lnLikelihood_rules ) );
        }

    };
}

#endif
