#include "RlTrace.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RlUserInterface.h"
#include "RlUtils.h"
#include "Integer.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "TypeSpec.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

Trace::Trace() : WorkspaceToCoreWrapperObject<RevBayesCore::TraceNumeric>()
{

    initMethods();

}


Trace::Trace(const RevBayesCore::TraceNumeric &t) : WorkspaceToCoreWrapperObject<RevBayesCore::TraceNumeric>( new RevBayesCore::TraceNumeric( t ) )
{

    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Trace* Trace::clone(void) const
{
    
	return new Trace(*this);
}


void Trace::computeStatistics( void )
{
    value->computeStatistics();
}


void Trace::constructInternalObject( void )
{
    
    throw RbException("We do not support a constructor function for Trace.");
}


/* Map calls to member methods */
RevPtr<RevVariable> Trace::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    if ( name == "setBurnin" )
    {
        found = true;

        int burnin = 0;

        RevObject& b = args[0].getVariable()->getRevObject();
        if ( b.isType( Integer::getClassTypeSpec() ) )
        {
            burnin = (int)static_cast<const Integer &>(b).getValue();
        }
        else
        {
            double burninFrac = static_cast<const Probability &>(b).getValue();
            burnin = int( floor( value->size()*burninFrac ) );
        }

        this->value->setBurnin( burnin );

        return NULL;
    }
    else if ( name == "size" || name == "getNumberSamples" )
    {
        found = true;

        bool post = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();

        int n = (int)this->value->size(post);

        return new RevVariable( new Natural( n ) );
    }
    else if ( name == "getBurnin" )
    {
        found = true;

        int n = this->value->getBurnin();

        return new RevVariable( new Natural( n ) );
    }
    else if (name == "summarize")
    {
        found = true;

        std::stringstream ss;
        ss << "Trace (";
        ss << getValue().getParameterName();
        ss << ")" << std::endl;
        ss << "Filename:                             " << value->getFileName() << std::endl;
        ss << "Parameter:                            " << value->getParameterName() << std::endl;
        ss << "Samples:                              " << value->size() << std::endl;
        ss << "ESS:                                  " << value->getESS() << std::endl;
        ss << "Burnin:                               " << value->getBurnin() << std::endl;
        ss << "Mean:                                 " << value->getMean() << std::endl;
        ss << "SEM:                                  " << value->getSEM() << std::endl;
        ss << "ESS > 625:                            " << (value->hasPassedEssThreshold() ? "TRUE" : "FALSE") << std::endl;
        ss << "PSRF < 1.01 (n = 10):                 " << (value->hasPassedGelmanRubinTest() ? "TRUE" : "FALSE") << std::endl;
        ss << "Geweke test (p = 0.01):               " << (value->hasPassedGewekeTest() ? "TRUE" : "FALSE") << std::endl;
        ss << "Stationarity test (p = 0.01, n = 10): " << (value->hasPassedStationarityTest() ? "TRUE" : "FALSE") << std::endl;

        RBOUT(ss.str());

        return NULL;
    }
    else if ( name == "getValues" )
    {
        found = true;
        
        const std::vector<double> &vals_tmp = value->getValues();
        size_t burnin = value->getBurnin();
        std::vector<double> vals = std::vector<double>(vals_tmp.size() - burnin, 0);
        for (size_t i=burnin; i<vals_tmp.size(); ++i)
        {
            vals[i-burnin] = vals_tmp[i];
        }
        bool positive = true;
        for (size_t i=0; i<vals.size(); ++i)
        {
            if ( vals[i] < 0.0 )
            {
                positive = false;
                break;
            }
        }
        if ( positive == true )
        {
            ModelVector<RealPos> *rl_vals = new ModelVector<RealPos>( vals );
            return new RevVariable( rl_vals );
        }
        else
        {
            ModelVector<Real> *rl_vals = new ModelVector<Real>( vals );
            return new RevVariable( rl_vals );
        }
        
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& Trace::getClassType(void)
{
    
    static std::string rev_type = "Trace";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Trace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::TraceNumeric>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}



/** Return member rules (no members) */
const MemberRules& Trace::getParameterRules(void) const
{
    
    static MemberRules modelMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules_set = true;
    }
    
    return modelMemberRules;
}


/** Get type spec */
const TypeSpec& Trace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void Trace::initMethods( void )
{
    ArgumentRules* burninArgRules = new ArgumentRules();
    std::vector<TypeSpec> burninTypes;
    burninTypes.push_back( Probability::getClassTypeSpec() );
    burninTypes.push_back( Integer::getClassTypeSpec() );
    burninArgRules->push_back( new ArgumentRule("burnin",      burninTypes, "The fraction/number of samples to disregard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "setBurnin", RlUtils::Void, burninArgRules) );

    ArgumentRules* getBurninArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "getBurnin", Natural::getClassTypeSpec(), getBurninArgRules) );

    ArgumentRules* summarizeArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarizeArgRules) );

    ArgumentRules* get_values_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "getValues", RlUtils::Void, get_values_arg_rules) );

    ArgumentRules* getNumberSamplesArgRules = new ArgumentRules();
    getNumberSamplesArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "getNumberSamples", Natural::getClassTypeSpec(), getNumberSamplesArgRules) );

    ArgumentRules* getSizeArgRules = new ArgumentRules();
    getSizeArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), getSizeArgRules) );
}


/** Get type spec */
void Trace::printValue(std::ostream &o, bool user) const
{
    o << "Trace (";
    o << getValue().getParameterName();
    o << ")\n";
}


/** Set a member variable */
void Trace::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "xxx") {
        
    } 
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
