#include "RlModelTrace.h"

#include <string>

#include "ArgumentRules.h"
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

ModelTrace::ModelTrace() : WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>()
{
    initMethods();
}


ModelTrace::ModelTrace(const RevBayesCore::ModelTrace &t) : WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>( new RevBayesCore::ModelTrace( t ) )
{
    initMethods();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
ModelTrace* ModelTrace::clone(void) const
{
    
    return new ModelTrace(*this);
}


void ModelTrace::constructInternalObject( void )
{
    throw RbException("We do not support a constructor function for ModelTrace.");
}


/* Map calls to member methods */
RevPtr<RevVariable> ModelTrace::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
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
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& ModelTrace::getClassType(void) {
    
    static std::string rev_type = "ModelTrace";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& ModelTrace::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules (no members) */
const MemberRules& ModelTrace::getParameterRules(void) const {
    
    static MemberRules modelMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules_set = true;
    }
    
    return modelMemberRules;
}


/** Get type spec */
const TypeSpec& ModelTrace::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void ModelTrace::initMethods( void )
{
    ArgumentRules* burninArgRules = new ArgumentRules();
    std::vector<TypeSpec> burninTypes;
    burninTypes.push_back( Probability::getClassTypeSpec() );
    burninTypes.push_back( Integer::getClassTypeSpec() );
    burninArgRules->push_back( new ArgumentRule("burnin",      burninTypes, "The fraction/number of samples to disregard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "setBurnin", RlUtils::Void, burninArgRules) );

    ArgumentRules* getBurninArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "getBurnin", Natural::getClassTypeSpec(), getBurninArgRules) );

    ArgumentRules* getNumberSamplesArgRules = new ArgumentRules();
    getNumberSamplesArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "getNumberSamples", Natural::getClassTypeSpec(), getNumberSamplesArgRules) );

    ArgumentRules* getSizeArgRules = new ArgumentRules();
    getSizeArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), getSizeArgRules) );
}


/** Get type spec */
void ModelTrace::printValue(std::ostream &o) const {
    
    o << "Filename:                  " << value->getFileName() << std::endl;
    o << "Parameter:                 " << value->getParameterName() << std::endl;
    o << "#Samples:                  " << value->size() << std::endl;
}


/** Set a member variable */
void ModelTrace::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "xxx") {
        
    }
    else {
        RevObject::setConstParameter(name, var);
    }
}


std::ostream& RevLanguage::operator<<(std::ostream& o, const RevLanguage::ModelTrace& x) {
    //    o << x.getParameterName();
    o << "ModelTrace (";
    //    const std::vector<double>& values = x.getValues();
    //    for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it) {
    //        if ( it != values.begin() ) {
    //            o << ", ";
    //        }
    //        o << *it;
    //    }
    o << ")";
    
    return o;
}
