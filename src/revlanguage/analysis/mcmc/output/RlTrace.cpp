#include "RlTrace.h"

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "Natural.h"
#include "Probability.h"
#include "RlUserInterface.h"
#include "RlUtils.h"

using namespace RevLanguage;

Trace::Trace() : WorkspaceToCoreWrapperObject<RevBayesCore::TraceNumeric>()
{

    ArgumentRules* summarizeArgRules = new ArgumentRules();
    std::vector<TypeSpec> burninTypes;
    burninTypes.push_back( Probability::getClassTypeSpec() );
    burninTypes.push_back( Integer::getClassTypeSpec() );
    summarizeArgRules->push_back( new ArgumentRule("burnin", burninTypes, "The fraction/number of samples to discregard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25)) );
    methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarizeArgRules) );

}


Trace::Trace(const RevBayesCore::TraceNumeric &t) : WorkspaceToCoreWrapperObject<RevBayesCore::TraceNumeric>( new RevBayesCore::TraceNumeric( t ) )
{

    ArgumentRules* summarizeArgRules = new ArgumentRules();
    std::vector<TypeSpec> burninTypes;
    burninTypes.push_back( Probability::getClassTypeSpec() );
    burninTypes.push_back( Integer::getClassTypeSpec() );
    summarizeArgRules->push_back( new ArgumentRule("burnin", burninTypes, "The fraction/number of samples to discregard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25)) );
    methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarizeArgRules) );

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
    
    if (name == "summarize")
    {
        found = true;

        std::stringstream ss;
        ss << "Trace (";
        ss << getValue().getParameterName();
        ss << ")" << std::endl;
        ss << "Filename:                  " << value->getFileName() << std::endl;
        ss << "Parameter:                 " << value->getParameterName() << std::endl;
        ss << "Samples:                   " << value->size() << std::endl;
        ss << "ESS:                       " << value->getESS() << std::endl;
        ss << "Burnin:                    " << value->getBurnin() << std::endl;
        ss << "Mean:                      " << value->getMean() << std::endl;
        ss << "SEM:                       " << value->getSEM() << std::endl;
        ss << "Geweke test:               " << (value->hasPassedGewekeTest() ? "TRUE" : "FALSE") << std::endl;
        ss << "Stationarity test:         " << (value->hasPassedStationarityTest() ? "TRUE" : "FALSE") << std::endl;

        RBOUT(ss.str());

        return NULL;
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
