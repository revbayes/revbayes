#ifndef Transform_Mul_H
#define Transform_Mul_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "TransformedDistribution.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

namespace RevLanguage {

    template <class T>
    class Transform_Mul : public TypedDistribution<T> {

    public:
        Transform_Mul( void );
        virtual ~Transform_Mul();

        // Basic utility functions
        Transform_Mul*                                  clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TransformedDistribution*          createDistribution(void) const;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:
        RevPtr<const RevVariable>                       base_distribution;
	RevPtr<const RevVariable>                       lambda;
    };


    template <class T>
    Transform_Mul<T>::Transform_Mul() : TypedDistribution< T >(),
                                     base_distribution( NULL )
    {
        this->markAsTransform();
    }

    template <class T>
    Transform_Mul<T>::~Transform_Mul()
    {
    }

    template <class T>
    Transform_Mul<T>* Transform_Mul<T>::clone( void ) const
    {
        return new Transform_Mul(*this);
    }

    template <class T>
    RevBayesCore::TransformedDistribution* Transform_Mul<T>::createDistribution( void ) const
    {
        using namespace Transforms;
        namespace Core = RevBayesCore;

        // get the parameters
        const Distribution& rl_vp               = static_cast<const Distribution &>( base_distribution->getRevObject() );
        Core::TypedDistribution<double>* vp     = static_cast<Core::TypedDistribution<double>* >( rl_vp.createDistribution() );

        Core::TypedDagNode<double>* l           = static_cast<const T &>( lambda->getRevObject() ).getDagNode();

        Core::TransformedDistribution* dist = new Core::TransformedDistribution(*vp, mul_transform, mul_inverse, log_mul_prime, {l});

        delete vp;

        return dist;
    }


    /* Get Rev type of object */
    template <class T>
    const std::string& Transform_Mul<T>::getClassType(void)
    {

        static std::string rev_type = "Transform_Mul";

        return rev_type;
    }

    /* Get class type spec describing type of object */
    template <class T>
    const TypeSpec& Transform_Mul<T>::getClassTypeSpec(void)
    {
        static TypeSpec rev_type_spec ( TypedDistribution< T >::getClassTypeSpec() );

        return rev_type_spec;
    }


    /**
     * Get the Rev name for the distribution.
     * This name is used for the constructor and the distribution functions,
     * such as the density and random value function
     *
     * \return Rev name of constructor function.
     */
    template <class T>
    std::string Transform_Mul<T>::getDistributionFunctionName( void ) const
    {
        // create a distribution name variable that is the same for all instance of this class
        std::string d_name = "_mul";

        return d_name;
    }


    /** Return member rules (no members) */
    template <class T>
    const MemberRules& Transform_Mul<T>::getParameterRules(void) const
    {
        static MemberRules dist_member_rules;
        static bool rules_set = false;

        if ( rules_set == false )
        {
            std::vector<TypeSpec> distTypes;

            if (Real::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Real>::getClassTypeSpec());

            if (RealPos::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<RealPos>::getClassTypeSpec());

            if (Probability::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Probability>::getClassTypeSpec());

            dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
            dist_member_rules.push_back( new ArgumentRule( "lambda", T::getClassTypeSpec()   , "The amount muled to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

            rules_set = true;
        }

        return dist_member_rules;
    }


    template <class T>
    const TypeSpec& Transform_Mul<T>::getTypeSpec( void ) const
    {
        static TypeSpec ts = getClassTypeSpec();

        return ts;
    }



    /** Set a member variable */
    template <class T>
    void Transform_Mul<T>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
    {
        if ( name == "baseDistribution" )
        {
            base_distribution = var;
        }
        else if ( name == "lambda" )
        {
            lambda = var;
        }
        else
        {
            TypedDistribution< T >::setConstParameter(name, var);
        }
    }


}

#endif // Transform_Mul_H
