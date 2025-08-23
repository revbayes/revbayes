#ifndef Dist_PseudoData_h
#define Dist_PseudoData_h

#include "PseudoDist.h"
#include "RlTypedDistribution.h"
#include "RlPseudoData.h"
#include "RlPseudoObservation.h"

namespace RevLanguage
{

    template <typename T>
    class Dist_Pseudo: public TypedDistribution<PseudoObservation>
    {
    public:
        Dist_Pseudo() {}

        Dist_Pseudo<T>* clone() const
        {
            return new Dist_Pseudo<T>(*this);
        }

        //!< Get Rev type
        static const std::string&                       getClassType(void)
        {
            static std::string rev_type = "PseudoDist";

            return rev_type;
        }

        //!< Get class type spec
        static const TypeSpec&                          getClassTypeSpec(void)
        {
            static TypeSpec rev_type_spec ( TypedDistribution<PseudoObservation>::getClassTypeSpec() );

            return rev_type_spec;
        }

        //!< Get the type spec of the instance
        const TypeSpec&                                 getTypeSpec(void) const
        {
            static TypeSpec ts = getClassTypeSpec();

            return ts;
        }

        /**
         * Get the Rev name for the distribution.
         * This name is used for the constructor and the distribution functions,
         * such as the density and random value function
         *
         * \return Rev name of constructor function.
         */
        //!< Get the Rev-name for this distribution.
        std::string                                     getDistributionFunctionName(void) const
        {
            return "Pseudo";
        }

        const MemberRules&                              getParameterRules(void) const;


        // Distribution functions you have to override
        RevBayesCore::PseudoDist<typename T::valueType>*                    createDistribution() const override;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    private:
        RevPtr<const RevVariable> parameter;
        RevPtr<const RevVariable> pseudoData;
    };



    template <typename T>
    RevBayesCore::PseudoDist<typename T::valueType>* Dist_Pseudo<T>::createDistribution( void ) const
    {
        // get the parameters
        auto param  = static_cast<const T&>( parameter->getRevObject() ).getDagNode();
        auto pd     = static_cast<const PseudoData<T>&>( pseudoData->getRevObject() ).getDagNode();

        return new RevBayesCore::PseudoDist<typename T::valueType>(param, pd);
    }

    template <typename T>
    const MemberRules&   Dist_Pseudo<T>::getParameterRules(void) const
    {
        static MemberRules dist_member_rules;
        static bool rules_set = false;

        if ( rules_set == false )
        {
            /* This should be simple but is a bit complicated because
             * the subtype relationships for T don't work on Distribution<T>
             *
             * That is:  RevBayes knows that RealPos < Real.
             *           RevBayes does NOT know that Distribution__RealPos < Distribution__Real.
             *
             * Therefore, depending on T, we need to list the types of the distributions that are allowed.
             * - for T=Real, the distribution can be on Real, RealPos, or Probability.
             * - for T=RealPos, the distribution can be on RealPos or Probability.
             * This encodes the fact that a RealPos plus a RealPos yields a RealPos.
             */
             
            dist_member_rules.push_back( new ArgumentRule( "parameter", T::getClassTypeSpec(), "The parameter of the unspecified distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
            dist_member_rules.push_back( new ArgumentRule( "pseudoData", PseudoData<T>::getClassTypeSpec(), "The parameter of the unspecified distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

            rules_set = true;
        }

        return dist_member_rules;
    }


    template <typename T>
    void Dist_Pseudo<T>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
    {
        if ( name == "parameter" )
        { 
            parameter = var;
        }
        else if ( name == "pseudoData" )
        { 
            pseudoData = var;
        }
        else
        {
            TypedDistribution< PseudoObservation >::setConstParameter(name, var);
        }
    }



}


#endif
