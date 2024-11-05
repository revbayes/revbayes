#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DistributionBeta.h"
#include "Mcmc.h"
#include "ModelVector.h"
#include "Natural.h"
#include "PowerPosteriorAnalysis.h"
#include "Probability.h"
#include "RevObject.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlPowerPosteriorAnalysis.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Monitor.h"
#include "Move.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevBayesCore { class Model; }


using namespace RevLanguage;

PowerPosteriorAnalysis::PowerPosteriorAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::PowerPosteriorAnalysis>()
{

    ArgumentRules* run_arg_rules = new ArgumentRules();
    run_arg_rules->push_back( new ArgumentRule("generations", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    run_arg_rules->push_back( new ArgumentRule("burninFraction", Probability::getClassTypeSpec(), "The fraction of samples to discard.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
    run_arg_rules->push_back( new ArgumentRule("preburninGenerations", Natural::getClassTypeSpec(), "The number of generations to run as pre-burnin when parameter tuning is done.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
    run_arg_rules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, run_arg_rules) );
    
    ArgumentRules* run_one_stone_arg_rules = new ArgumentRules();
    run_one_stone_arg_rules->push_back( new ArgumentRule("index", Natural::getClassTypeSpec(), "Index of the stone/power to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    run_one_stone_arg_rules->push_back( new ArgumentRule("generations", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    run_one_stone_arg_rules->push_back( new ArgumentRule("burninFraction", Probability::getClassTypeSpec(), "The fraction of samples to discard.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
    run_one_stone_arg_rules->push_back( new ArgumentRule("preburninGenerations", Natural::getClassTypeSpec(), "The number of generations to run as pre-burnin when parameter tuning is done.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
    run_one_stone_arg_rules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
    methods.addFunction( new MemberProcedure( "runOneStone", RlUtils::Void, run_one_stone_arg_rules) );
    
    ArgumentRules* summarize_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarize_arg_rules) );

    ArgumentRules* burnin_arg_rules = new ArgumentRules();
    burnin_arg_rules->push_back( new ArgumentRule("generations"   , Natural::getClassTypeSpec(), "The number of generations to run.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    burnin_arg_rules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The frequency when the moves are tuned (usually between 50 and 1000).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "burnin", RlUtils::Void, burnin_arg_rules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
PowerPosteriorAnalysis* PowerPosteriorAnalysis::clone(void) const
{

	return new PowerPosteriorAnalysis(*this);
}


void PowerPosteriorAnalysis::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new sliding move
    const RevBayesCore::Model&                      mdl     = static_cast<const Model &>( model->getRevObject() ).getValue();
    const WorkspaceVector<Move>&                    rlmvs   = static_cast<const WorkspaceVector<Move> &>( moves->getRevObject() );
    const WorkspaceVector<Monitor>&                 rlmntr  = static_cast<const WorkspaceVector<Monitor> &>( monitors->getRevObject() );
    RevBayesCore::RbVector<RevBayesCore::Monitor>   mntr;
    for ( size_t i = 0; i < rlmntr.size(); ++i )
    {
        mntr.push_back( rlmntr[i].getValue() );
    }
    RevBayesCore::RbVector<RevBayesCore::Move>      mvs;
    for ( size_t i = 0; i < rlmvs.size(); ++i )
    {
        mvs.push_back( rlmvs[i].getValue() );
    }
    const std::string&                              fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const double                                    alpha   = static_cast<const RealPos &>( alphaVal->getRevObject() ).getValue();
    const int                                       sf      = (int)static_cast<const Natural &>( sampFreq->getRevObject() ).getValue();
    const int                                       k       = (int)static_cast<const Natural &>( proc_per_lik->getRevObject() ).getValue();

    RevBayesCore::Mcmc *m = new RevBayesCore::Mcmc(mdl, mvs, mntr);
    m->setScheduleType( "random" );

    value = new RevBayesCore::PowerPosteriorAnalysis( m, fn, size_t(k) );

    std::vector<double> beta;
    if ( powers->getRevObject() != RevNullObject::getInstance() )
    {
        beta = static_cast<const ModelVector<RealPos> &>( powers->getRevObject() ).getValue();
    }
    else
    {
        int k = (int)static_cast<const Natural &>( cats->getRevObject() ).getValue();
        for (int i = k; i >= 0; --i)
        {
            double b = RevBayesCore::RbStatistics::Beta::quantile(alpha,1.0,i / double(k));
            beta.push_back( b );
        }
    }

    value->setPowers( beta );
    value->setSampleFreq( sf );
}


/* Map calls to member methods */
RevPtr<RevVariable> PowerPosteriorAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{

    if (name == "run")
    {
        found = true;

        // get the member with give index
        std::int64_t gen = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        double burn_frac = static_cast<const Probability &>( args[1].getVariable()->getRevObject() ).getValue();
        size_t preburn_gen = gen;
        if ( args[2].getVariable()->getRevObject() != RevNullObject::getInstance() )
        {
            preburn_gen = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();
        }
        size_t tune_int = static_cast<const Natural &>( args[3].getVariable()->getRevObject() ).getValue();
        value->runAll( size_t(gen), burn_frac, preburn_gen, tune_int );

        return NULL;
    }
    else if (name == "runOneStone")
    {
        found = true;
        
        /* We will keep the C++ indexing from 0. This is un-Rev-like, but it is in keeping with the
         * 'cats' argument, which functions in the same way (i.e., specifying cats=127 will run 128
         * stones with indices ranging from 0 to 127)
         */
        size_t ind = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        std::int64_t gen = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
        double burn_frac = static_cast<const Probability &>( args[2].getVariable()->getRevObject() ).getValue();
        size_t preburn_gen = gen;
        if ( args[3].getVariable()->getRevObject() != RevNullObject::getInstance() )
        {
            preburn_gen = static_cast<const Natural &>( args[3].getVariable()->getRevObject() ).getValue();
        }
        size_t tune_int = static_cast<const Natural &>( args[4].getVariable()->getRevObject() ).getValue();
        value->runStone( ind, size_t(gen), burn_frac, preburn_gen, tune_int );
        
        return NULL;
    }
    else if (name == "summarize")
    {
        found = true;
        
        value->summarizeStones();
        
        return NULL;
    }
    else if (name == "burnin")
    {
        found = true;

        // get the member with give index
        int gen = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        int tuningInterval = (int)static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
        value->burnin( size_t(gen), size_t(tuningInterval) );

        return NULL;
    }

    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& PowerPosteriorAnalysis::getClassType(void)
{

    static std::string rev_type = "PowerPosteriorAnalysis";

	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& PowerPosteriorAnalysis::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::PowerPosteriorAnalysis>::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string PowerPosteriorAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "powerPosterior";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& PowerPosteriorAnalysis::getParameterRules(void) const
{

    static MemberRules member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {

        member_rules.push_back( new ArgumentRule("model"      , Model::getClassTypeSpec()                   , "The model graph.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("moves"      , WorkspaceVector<Move>::getClassTypeSpec()   , "The vector moves to use.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("monitors"   , WorkspaceVector<Monitor>::getClassTypeSpec(), "The monitors to call. Do not provide a screen monitor.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("filename"   , RlString::getClassTypeSpec()                , "The name of the file for the likelihood samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule("powers"     , ModelVector<RealPos>::getClassTypeSpec()    , "A vector of powers.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        member_rules.push_back( new ArgumentRule("cats"       , Natural::getClassTypeSpec()                 , "The number of categories if no powers are specified.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
        member_rules.push_back( new ArgumentRule("alpha"      , RealPos::getClassTypeSpec()                 , "The alpha parameter of the beta distribution if no powers are specified.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.2) ) );
        member_rules.push_back( new ArgumentRule("sampleFreq" , Natural::getClassTypeSpec()                 , "The sampling frequency of the likelihood values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
        member_rules.push_back( new ArgumentRule("procPerLikelihood" , Natural::getClassTypeSpec()          , "Number of processors used to compute the likelihood.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        rules_set = true;
    }

    return member_rules;
}


/** Get type spec */
const TypeSpec& PowerPosteriorAnalysis::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void PowerPosteriorAnalysis::printValue(std::ostream &o) const
{

    o << "PowerPosterior";
}


/** Set a member variable */
void PowerPosteriorAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "model")
    {
        model = var;
    }
    else if ( name == "moves")
    {
        moves = var;
    }
    else if ( name == "monitors")
    {
        monitors = var;
    }
    else if ( name == "filename")
    {
        filename = var;
    }
    else if ( name == "cats")
    {
        cats = var;
    }
    else if ( name == "powers")
    {
        powers = var;
    }
    else if ( name == "alpha")
    {
        alphaVal = var;
    }
    else if ( name == "sampleFreq")
    {
        sampFreq = var;
    }
    else if ( name == "procPerLikelihood" )
    {
        proc_per_lik = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}
