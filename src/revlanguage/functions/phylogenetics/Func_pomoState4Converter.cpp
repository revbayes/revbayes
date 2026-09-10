#include "Func_pomoState4Converter.h"

#include <map>
#include <string>
#include <vector>

#include "ModelVector.h"
#include "Natural.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTaxon.h"
#include "PoMoState4Converter.h"
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
Func_pomoState4Converter::Func_pomoState4Converter( void ) : Procedure( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_pomoState4Converter* Func_pomoState4Converter::clone( void ) const {
    
    return new Func_pomoState4Converter( *this );
}


RevPtr<RevVariable> Func_pomoState4Converter::execute() {
    
    const RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* aln = static_cast<const AbstractHomologousDiscreteCharacterData&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode< std::int64_t >* n = static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getDagNode();


    RevBayesCore::PoMoState4Converter* c = new RevBayesCore::PoMoState4Converter(  );
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::Taxon> >* taxa  = static_cast< const ModelVector<Taxon> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    std::map <std::string, std::string > gene2species;

    for (RevBayesCore::RbIterator<RevBayesCore::Taxon> it=taxa->getValue().begin(); it!=taxa->getValue().end(); ++it)
    {
        gene2species[it->getName()] = it->getSpeciesName();
    }
    
    AbstractHomologousDiscreteCharacterData PoMoAln = c->convertData( aln->getValue(), (int)n->getValue(), gene2species ) ;
        
    return new RevVariable( new AbstractHomologousDiscreteCharacterData( PoMoAln ) );
}


/* Get argument rules */
const ArgumentRules& Func_pomoState4Converter::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "aln"      , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "virtualNe", Natural::getClassTypeSpec()                                , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec()                     , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_pomoState4Converter::getClassType(void)
{
    
    static std::string rev_type = "Func_pomoState4Converter";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_pomoState4Converter::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_pomoState4Converter::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "pomoState4Convert";
    
    return f_name;
}


/* Get return type */
const TypeSpec& Func_pomoState4Converter::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = AbstractHomologousDiscreteCharacterData::getClassTypeSpec();
    
    return return_typeSpec;
}


const TypeSpec& Func_pomoState4Converter::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
