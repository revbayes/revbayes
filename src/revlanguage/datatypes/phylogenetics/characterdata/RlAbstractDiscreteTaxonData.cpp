#include "RlAbstractDiscreteTaxonData.h"

#include <cstddef>

#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RlDiscreteCharacterState.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "MethodTable.h"
#include "RbException.h"
#include "RevVariable.h"
#include "TypeSpec.h"

namespace RevLanguage { class RevObject; }


AbstractDiscreteTaxonData::AbstractDiscreteTaxonData(void) : ModelObject<RevBayesCore::AbstractDiscreteTaxonData >()
{
    
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", DiscreteCharacterState::getClassTypeSpec(), squareBracketArgRules) );
    
    // Add method for call "nchar" as a function
    ArgumentRules* nchar_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "nchar", Natural::getClassTypeSpec(), nchar_arg_rules) );
    
}


AbstractDiscreteTaxonData::AbstractDiscreteTaxonData( RevBayesCore::AbstractDiscreteTaxonData *v) : ModelObject<RevBayesCore::AbstractDiscreteTaxonData >( v )
{

    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", DiscreteCharacterState::getClassTypeSpec(), squareBracketArgRules) );

    // Add method for call "nchar" as a function
    ArgumentRules* nchar_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "nchar", Natural::getClassTypeSpec(), nchar_arg_rules) );
    
}


AbstractDiscreteTaxonData::AbstractDiscreteTaxonData( const RevBayesCore::AbstractDiscreteTaxonData &d) : ModelObject<RevBayesCore::AbstractDiscreteTaxonData >( d.clone() )
{

    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", DiscreteCharacterState::getClassTypeSpec(), squareBracketArgRules) );

    // Add method for call "nchar" as a function
    ArgumentRules* nchar_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "nchar", Natural::getClassTypeSpec(), nchar_arg_rules) );
    
}


AbstractDiscreteTaxonData::AbstractDiscreteTaxonData( RevBayesCore::TypedDagNode<RevBayesCore::AbstractDiscreteTaxonData> *d) : ModelObject<RevBayesCore::AbstractDiscreteTaxonData >( d )
{

    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", DiscreteCharacterState::getClassTypeSpec(), squareBracketArgRules) );

    // Add method for call "nchar" as a function
    ArgumentRules* nchar_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "nchar", Natural::getClassTypeSpec(), nchar_arg_rules) );
    
}



AbstractDiscreteTaxonData* AbstractDiscreteTaxonData::clone() const
{
    return new AbstractDiscreteTaxonData( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> AbstractDiscreteTaxonData::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if ( name == "[]")
    {
        found = true;
        
        // get the member with give index
        const Natural &index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() );
        
        if (this->dag_node->getValue().getNumberOfCharacters() < (size_t)(index.getValue()) )
        {
            throw RbException("Index out of bounds in []");
        }
        
        RevObject* element = new DiscreteCharacterState( this->dag_node->getValue().getCharacter( size_t(index.getValue()) - 1) );
        return new RevVariable( element );
    }
    else if ( name == "nchar")
    {
        found = true;
            
        // get the member with give index
        size_t n_char = this->dag_node->getValue().getNumberOfCharacters();
            
        Natural* element = new Natural( n_char );
        return new RevVariable( element );
    }
    
    return ModelObject<RevBayesCore::AbstractDiscreteTaxonData >::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& AbstractDiscreteTaxonData::getClassType(void)
{
    
    static std::string rev_type = "AbstractDiscreteTaxonData";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& AbstractDiscreteTaxonData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::AbstractDiscreteTaxonData >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& AbstractDiscreteTaxonData::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
