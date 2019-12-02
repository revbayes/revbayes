#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DistributionBeta.h"
#include "Mcmc.h"
#include "ModelVector.h"
#include "Model.h"
#include "Natural.h"
#include "AdaptivePowerPosteriorAnalysis.h"
#include "RevObject.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlModel.h"
#include "RlMonitor.h"
#include "RlMove.h"
#include "RlAdaptivePowerPosteriorAnalysis.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"


using namespace RevLanguage;

AdaptivePowerPosteriorAnalysis::AdaptivePowerPosteriorAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::AdaptivePowerPosteriorAnalysis>()
{

    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule("generations",    Natural::getClassTypeSpec(), "The number of (post-burnin) generations to run for each stone.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    runArgRules->push_back( new ArgumentRule("burnin",         Natural::getClassTypeSpec(), "The number of burnin generations to run for each stone.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    runArgRules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The frequency when the moves are tuned (usually between 50 and 1000) during the burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );

    ArgumentRules* burninArgRules = new ArgumentRules();
    burninArgRules->push_back( new ArgumentRule("generations"   , Natural::getClassTypeSpec(), "The number of generations to run for each stone.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    burninArgRules->push_back( new ArgumentRule("tuningInterval", Natural::getClassTypeSpec(), "The frequency when the moves are tuned (usually between 50 and 1000) during the burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "burnin", RlUtils::Void, burninArgRules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
AdaptivePowerPosteriorAnalysis* AdaptivePowerPosteriorAnalysis::clone(void) const
{

	return new AdaptivePowerPosteriorAnalysis(*this);
}


void AdaptivePowerPosteriorAnalysis::constructInternalObject( void )
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
    const int                                       sf      = (int)static_cast<const Natural &>( sampFreq->getRevObject() ).getValue();
    const double                                    tol     = (double)static_cast<const RealPos &>( tolerance->getRevObject() ).getValue();
    const int                                       k       = (int)static_cast<const Natural &>( proc_per_lik->getRevObject() ).getValue();

    RevBayesCore::Mcmc *m = new RevBayesCore::Mcmc(mdl, mvs, mntr);
    m->setScheduleType( "random" );

    value = new RevBayesCore::AdaptivePowerPosteriorAnalysis( m, fn, size_t(k) );

    value->setTolerance( tol );
    value->setSampleFreq( sf );
    
}


/* Map calls to member methods */
RevPtr<RevVariable> AdaptivePowerPosteriorAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{

    if (name == "run")
    {

        found = true;

        // get the member with give index
        int gen  = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        int burn = (int)static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
        int ti   = (int)static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue();

        value->runAll( size_t(gen), size_t(burn), size_t(ti) );

        
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
const std::string& AdaptivePowerPosteriorAnalysis::getClassType(void)
{

    static std::string rev_type = "AdaptivePowerPosteriorAnalysis";

	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& AdaptivePowerPosteriorAnalysis::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::AdaptivePowerPosteriorAnalysis>::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string AdaptivePowerPosteriorAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "adaptivePowerPosterior";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& AdaptivePowerPosteriorAnalysis::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        memberRules.push_back( new ArgumentRule("model"      , Model::getClassTypeSpec()                   , "The model graph.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("moves"      , WorkspaceVector<Move>::getClassTypeSpec()   , "The vector moves to use.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("monitors"   , WorkspaceVector<Monitor>::getClassTypeSpec(), "The monitors to call. Do not provide a screen monitor.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("filename"   , RlString::getClassTypeSpec()                , "The name of the file for the likelihood samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("sampleFreq" , Natural::getClassTypeSpec()                 , "The sampling frequency of the likelihood values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(100) ) );
        memberRules.push_back( new ArgumentRule("tol"        , RealPos::getClassTypeSpec()                 , "The tolerance for stopping the adaptive sampler.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.1) ) );
        memberRules.push_back( new ArgumentRule("procPerLikelihood" , Natural::getClassTypeSpec()          , "Number of processors used to compute the likelihood.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        rules_set = true;
    }

    return memberRules;
}


/** Get type spec */
const TypeSpec& AdaptivePowerPosteriorAnalysis::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void AdaptivePowerPosteriorAnalysis::printValue(std::ostream &o) const
{

    o << "AdaptivePowerPosterior";
}


/** Set a member variable */
void AdaptivePowerPosteriorAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "sampleFreq")
    {
        sampFreq = var;
    }
    else if ( name == "tol")
    {
        tolerance = var;
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
