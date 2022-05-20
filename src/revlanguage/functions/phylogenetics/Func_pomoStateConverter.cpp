#include "Func_pomoStateConverter.h"

#include <map>
#include <string>
#include <vector>

#include "ModelVector.h"
#include "Natural.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTaxon.h"
#include "PoMoStateConverter.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "HomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TypeSpec.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }

using namespace RevLanguage;

/** default constructor */
Func_pomoStateConverter::Func_pomoStateConverter( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_pomoStateConverter* Func_pomoStateConverter::clone( void ) const
{
    
    return new Func_pomoStateConverter( *this );
}


RevPtr<RevVariable> Func_pomoStateConverter::execute()
{
    
    const RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* aln = static_cast<const AbstractHomologousDiscreteCharacterData&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    long num_states     = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();
    long virt_pop_size  = static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getValue();


    RevBayesCore::PoMoStateConverter* c = new RevBayesCore::PoMoStateConverter(  );
    
    const RevBayesCore::RbVector<RevBayesCore::Taxon>& taxa  = static_cast< const ModelVector<Taxon> &>( this->args[3].getVariable()->getRevObject() ).getValue();
    
    AbstractHomologousDiscreteCharacterData* PoMoAln = NULL;
    
    if ( num_states == 2 )
    {
        PoMoAln = new AbstractHomologousDiscreteCharacterData( c->convertData2( aln->getValue(), virt_pop_size, taxa ) );
    }
    else if ( num_states == 4 )
    {
//        c->convertData4(aln->getValue(), virt_pop_size, taxa);
    }
    return new RevVariable( PoMoAln );
}


/* Get argument rules */
const ArgumentRules& Func_pomoStateConverter::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argument_rules.push_back( new ArgumentRule( "aln"      , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The traditional multiple sequence alignment of all sequences.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "k"        , Natural::getClassTypeSpec()                                , "The number of states.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "virtualNe", Natural::getClassTypeSpec()                                , "The virtual population size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec()                     , "The taxa to match the individuals to species/populations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argument_rules;
}


const std::string& Func_pomoStateConverter::getClassType(void)
{
    
    static std::string rev_type = "Func_pomoStateConverter";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_pomoStateConverter::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_pomoStateConverter::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pomoStateConvert";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_pomoStateConverter::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_pomoStateConverter::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
