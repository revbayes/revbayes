#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "Natural.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Constructor from empty RlString */
RlString::RlString( void ) : ModelObject<std::string>()
{
    
    initMethods();
}


/** Constructor from RlString */
RlString::RlString(const std::string& v, bool parse) : ModelObject<std::string>( new std::string(v) )
{
    
    if ( parse == true )
    {
        parseValue();
    }
    initMethods();
}


/* Construct from DAG node */
RlString::RlString( RevBayesCore::TypedDagNode<std::string> *v, bool parse ) : ModelObject<std::string>( v )
{
    
    if ( parse == true )
    {
        parseValue();
    }
    initMethods();
}


/**
 * Generic addition operator.
 * We test if the rhs is of a type that we use for a specialized addition operation.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RevObject* RlString::add( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec() == RlString::getClassTypeSpec() )
    {
        return add( static_cast<const RlString&>( rhs ) );
    }
    
    return ModelObject<std::string>::add( rhs );
}


/**
 * Specialized addition operation between two string values.
 * The return value is also of type string.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the concatenation.
 */
RlString* RlString::add(const RevLanguage::RlString &rhs) const
{
    
    RlString *n = new RlString( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/** Clone function */
RlString* RevLanguage::RlString::clone() const
{

	return new RlString(*this);
}


/**
 * Map calls to member methods.
 */
RevPtr<RevVariable> RlString::executeMethod( std::string const &name, const std::vector<Argument> &args, bool &found )
{
    
    if ( name == "size" )
    {
        found = true;
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new Natural( getValue().size() ), "" ) );
    }
    else if ( name == "charAt" )
    {
        found = true;
        
        long index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        const std::string &str = getValue();
        std::string letter(1,str[index]);
        return RevPtr<RevVariable>( new RevVariable( new RlString( letter ) ) );
    }
    else if ( name == "substr" )
    {
        found = true;
        
        size_t arg_idx = 0;
        long begin = static_cast<const Natural&>( args[arg_idx++].getVariable()->getRevObject() ).getValue() - 1;
        long end = static_cast<const Natural&>( args[arg_idx++].getVariable()->getRevObject() ).getValue() - 1;
        const std::string &str = getValue();
        std::string substr = str.substr(begin,end-begin+1);
        return RevPtr<RevVariable>( new RevVariable( new RlString( substr ) ) );
    }
    else if ( name == "indexOf" || name == "find" )
    {
        found = true;
        
        size_t arg_idx = 0;
        const std::string substr = static_cast<const RlString&>( args[arg_idx++].getVariable()->getRevObject() ).getValue();

        std::string val = getValue();
        size_t index = val.find( substr ) + 1;
        return RevPtr<RevVariable>( new RevVariable( new Natural( index ) ) );
    }
    else if ( name == "[]" )
    {
        found = true;
        
        long index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        const std::string &str = getValue();
        std::string letter(1,str[index]);
        return RevPtr<RevVariable>( new RevVariable( new RlString( letter ) ) );
    }
    
    return ModelObject<std::string>::executeMethod( name, args, found );
}



/** Get Rev type of object */
const std::string& RlString::getClassType(void)
{
    
    static std::string rev_type = "String";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& RlString::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& RlString::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
  * Initialize the methods.
  */
void RlString::initMethods( void )
{
    
    ArgumentRules* size_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), size_arg_rules) );

    ArgumentRules* char_at_arg_rules = new ArgumentRules();
    char_at_arg_rules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "charAt", RlString::getClassTypeSpec(), char_at_arg_rules) );

    ArgumentRules* substr_arg_rules = new ArgumentRules();
    substr_arg_rules->push_back( new ArgumentRule( "begin", Natural::getClassTypeSpec(), "The index of the first character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    substr_arg_rules->push_back( new ArgumentRule( "end",   Natural::getClassTypeSpec(), "The index of the last character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "substr", RlString::getClassTypeSpec(), substr_arg_rules) );

    ArgumentRules* index_of_arg_rules = new ArgumentRules();
    index_of_arg_rules->push_back( new ArgumentRule( "substr", RlString::getClassTypeSpec(), "The substring for which we want to find the position.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "indexOf", Natural::getClassTypeSpec(), index_of_arg_rules) );

    ArgumentRules* find_arg_rules = new ArgumentRules();
    find_arg_rules->push_back( new ArgumentRule( "substr", RlString::getClassTypeSpec(), "The substring for which we want to find the position.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "find", Natural::getClassTypeSpec(), find_arg_rules) );

    ArgumentRules* element_arg_rules = new ArgumentRules();
    element_arg_rules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", RlString::getClassTypeSpec(), element_arg_rules ) );
    
}



/** Print value */
void RlString::parseValue(void)
{
    
    const std::string &v = dag_node->getValue();
    
    std::string res;
    std::string::const_iterator it = v.begin();
    while (it != v.end())
    {
        char c = *it++;
        if (c == '\\' && it != v.end())
        {
            switch (*it++)
            {
                case '\\': c = '\\'; break;
                case 'n': c = '\n'; break;
                case 't': c = '\t'; break;
                case '"': c = '"'; break;
                    // all other escapes
                default:
                    // invalid escape sequence - skip it. alternatively you can copy it as is, throw an exception...
                    continue;
            }
        }
        res += c;
    }
    
    this->setValue( new std::string(res) );
}

