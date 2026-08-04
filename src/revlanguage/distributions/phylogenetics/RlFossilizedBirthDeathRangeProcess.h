#ifndef RlFossilizedBirthDeathRangeProcess_H
#define RlFossilizedBirthDeathRangeProcess_H

#include "DistributionMemberFunction.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDistributionMemberFunction.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RlDistribution.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the FossilizedBirthDeathRangeProcess
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    template <typename rlType>
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<rlType> {
        
    public:
        virtual                                             ~FossilizedBirthDeathRangeProcess(void);                                          //!< Destructor
        
        // Basic utility functions you have to overwrite
        virtual FossilizedBirthDeathRangeProcess<rlType>*   clone(void) const = 0;                                                              //!< Clone the object

        // Basic utility functions you may want to overwrite
        const MemberRules&                                  getParameterRules(void) const;                                                      //!< Get member rules (const), including the reporting args
        MethodTable                                         getDistributionMethods(void) const;                                                 //!< The range times and appearances, which no monitor can otherwise reach
        
        // Basic utility functions
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec

    protected:
        FossilizedBirthDeathRangeProcess( void );

        const MemberRules&                                  getCoreParameterRules(void) const;                                              //!< Get member rules (const), without the reporting args
        static void                                         appendParameterRules(MemberRules &rules, bool include_reporting);                   //!< Build the member rules, with or without the reporting args
        static ArgumentRule*                                originPriorRule(void);                                                              //!< The optional origin prior argument, for the processes that take one
        RevBayesCore::TypedDistribution<double>*            createOriginPrior(void) const;                                                      //!< Build the core origin prior, or NULL if none was given
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
    
        // members        
        RevPtr<const RevVariable>                           lambda;                                                                             //!< The speciation rate(s)
        RevPtr<const RevVariable>                           mu;                                                                                 //!< The extinction rate(s)
        RevPtr<const RevVariable>                           psi;                                                                                //!< The fossilization rate(s)
        RevPtr<const RevVariable>                           rho;                                                                                //!< The extant sampling proportion
        RevPtr<const RevVariable>                           timeline;                                                                           //!< The interval times
        RevPtr<const RevVariable>                           present;                                                                            //!< Where the process stops
        RevPtr<const RevVariable>                           taxa;                                                                               //!< The taxa
        RevPtr<const RevVariable>                           condition;                                                                          //!< The condition of the process
        RevPtr<const RevVariable>                           complete;                                                                           //!< Is the fossil record complete?
        RevPtr<const RevVariable>                           origin_prior;                                                                       //!< Optional prior on the origin, which is the oldest birth

    };
    
}


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::FossilizedBirthDeathRangeProcess() : TypedDistribution<rlType>()
{

}


/**
 * Default destructor.
 *
 * The default destructor does nothing.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::~FossilizedBirthDeathRangeProcess()
{

}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
template <typename rlType>
const std::string& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassType(void)
{

    static std::string rev_type = "FossilizedBirthDeathRange";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
template <typename rlType>
const TypeSpec& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<rlType>::getClassTypeSpec() ) );

    return rev_type_spec;
}



/**
 * Build the member rules shared by the matrix and tree range processes.
 *
 * The reporting args (complete, reporting) belong to the fossil-record term, so the bare range process
 * omits them and takes its reporting model from a dnFossilRecord node instead. Keep them at this
 * position in the list: the argument order is part of the interface.
 *
 * \param[out]   rules              The rule list to append to.
 * \param[in]    include_reporting  Whether to include the reporting args.
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::appendParameterRules(MemberRules &rules, bool include_reporting)
{
    std::vector<TypeSpec> paramTypes;
    paramTypes.push_back( RealPos::getClassTypeSpec() );
    paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
    rules.push_back( new ArgumentRule( "lambda",  paramTypes, "The speciation rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    rules.push_back( new ArgumentRule( "mu",      paramTypes, "The extinction rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "psi",     paramTypes, "The fossil sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "rho",     Probability::getClassTypeSpec(), "The extant sampling fraction.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

    rules.push_back( new ArgumentRule( "timeline",   ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
    rules.push_back( new ArgumentRule( "present",    RealPos::getClassTypeSpec(), "The time defining the present. Minimum age of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

    std::vector<std::string> optionsCondition;
    optionsCondition.push_back( "time" );
    optionsCondition.push_back( "sampling" );
    optionsCondition.push_back( "survival" );
    rules.push_back( new OptionRule( "condition", new RlString("time"), optionsCondition, "The condition of the process." ) );
    rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa with fossil occurrence information.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    if ( include_reporting )
    {
        rules.push_back( new ArgumentRule( "complete", RlBoolean::getClassTypeSpec(), "Is the fossil record complete (every sampled occurrence reported)? FALSE is the first/last rule: both extremes reported, the interior count marginalized.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
    }

}


/**
 * The optional origin prior argument. The origin is the oldest birth, so the prior is evaluated
 * there rather than on a parameter of its own.
 */
template <typename rlType>
ArgumentRule* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::originPriorRule(void)
{
    return new ArgumentRule( "origin_prior", TypedDistribution<RealPos>::getClassTypeSpec(), "A prior on the origin of the process, which is the oldest birth.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL );
}


template <typename rlType>
RevBayesCore::TypedDistribution<double>* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::createOriginPrior(void) const
{
    if ( origin_prior == NULL || origin_prior->getRevObject() == RevNullObject::getInstance() )
    {
        return NULL;
    }

    const Distribution &rl_op = static_cast<const Distribution &>( origin_prior->getRevObject() );

    return static_cast<RevBayesCore::TypedDistribution<double>* >( rl_op.createDistribution() );
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the fossilized birth-death range process are:
 * (1) time of the process since the origin.
 * (2) time of the process since the rootAge.
 * (3) the sampling probability.
 * (4) the sampling strategy.
 * (5) the condition.
 * (6) the number of taxa.
 * (7) the taxon names.
 * (8) the clade constraints.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules, true );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Get the member rules of the bare b/d process: everything except the reporting args.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getCoreParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules, false );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "psi" )
    {
        psi = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "present" )
    {
        present = var;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "complete" )
    {
        complete = var;
    }
    else if ( name == "origin_prior" )
    {
        origin_prior = var;
    }
    else
    {
        TypedDistribution<rlType>::setConstParameter(name,var);
    }

}


/** The range times and appearances, which are internal to the distribution. */
template <typename rlType>
RevLanguage::MethodTable RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<rlType>::getDistributionMethods();

    ArgumentRules* first_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getFirstAppearances", this->variable, first_app_arg_rules, true ) );

    ArgumentRules* last_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getLastAppearances", this->variable, last_app_arg_rules, true ) );

    ArgumentRules* origin_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, RealPos >( "getOrigin", this->variable, origin_arg_rules, true ) );

    ArgumentRules* origination_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getOriginationTimes", this->variable, origination_arg_rules, true ) );

    ArgumentRules* extinction_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getExtinctionTimes", this->variable, extinction_arg_rules, true ) );

    return methods;
}

#endif
