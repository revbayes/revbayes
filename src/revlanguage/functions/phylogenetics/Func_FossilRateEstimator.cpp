#include "Func_FossilRateEstimator.h"

#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "FossilIntervalCounts.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlMatrixReal.h"
#include "RlString.h"
#include "RevNullObject.h"
#include "RlTaxon.h"

using namespace RevLanguage;


Func_FossilRateEstimator::Func_FossilRateEstimator( void ) : TypedFunction< MatrixReal >()
{

}


const ArgumentRules& Func_FossilRateEstimator::getArgumentRules( void ) const
{
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {
        argument_rules.push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "The taxa with fossil occurrence information.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "timeline", ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        argument_rules.push_back( new ArgumentRule( "present", RealPos::getClassTypeSpec(), "The time defining the present. Minimum age of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        argument_rules.push_back( new ArgumentRule( "max_age", RealPos::getClassTypeSpec(), "The older edge of the oldest interval. Left out, that interval is unbounded and has no rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<std::string> ambiguous_options;
        ambiguous_options.push_back( "midpoint" );
        ambiguous_options.push_back( "overlap" );
        ambiguous_options.push_back( "exclude" );
        argument_rules.push_back( new OptionRule( "ambiguous", new RlString("midpoint"), ambiguous_options, "Where an occurrence whose reported bin straddles an interval boundary is counted. Taxa only." ) );

        rules_set = true;
    }

    return argument_rules;
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_FossilRateEstimator::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* timeline =
        static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode< double >* present =
        static_cast<const RealPos &>( this->args[2].getVariable()->getRevObject() ).getDagNode();

    const RevBayesCore::TypedDagNode< double >* max_age = NULL;
    if ( this->args[3].getVariable()->getRevObject() != RevNullObject::getInstance() )
    {
        max_age = static_cast<const RealPos &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    }

    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::Taxon> >* taxa =
        static_cast<const ModelVector<Taxon> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    const std::string &ambiguous = static_cast<const RlString &>( this->args[4].getVariable()->getRevObject() ).getValue();
    RevBayesCore::FossilAgeAmbiguity policy = RevBayesCore::FOSSIL_AGE_MIDPOINT;
    if ( ambiguous == "overlap" )      policy = RevBayesCore::FOSSIL_AGE_OVERLAP;
    else if ( ambiguous == "exclude" ) policy = RevBayesCore::FOSSIL_AGE_EXCLUDE;

    return new RevBayesCore::FossilRateEstimatorFunction( taxa, timeline, present, max_age, method(), policy );
}
