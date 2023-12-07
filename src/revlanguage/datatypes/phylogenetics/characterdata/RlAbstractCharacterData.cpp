#include "RlAbstractCharacterData.h"

#include <cstddef>
#include <iostream>
#include <string>

#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "AbstractCharacterData.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

AbstractCharacterData::AbstractCharacterData( void )
{

}


MethodTable AbstractCharacterData::getCharacterDataMethods( void ) const
{
    
    MethodTable methods = MethodTable();
    
    ArgumentRules* addTaxonArgRules             = new ArgumentRules();
    ArgumentRules* namesArgRules                = new ArgumentRules();
    ArgumentRules* excludeTaxaArgRules          = new ArgumentRules();
    ArgumentRules* excludeTaxaArgRules2         = new ArgumentRules();
    ArgumentRules* filenameArgRules             = new ArgumentRules();
    ArgumentRules* includeTaxaArgRules          = new ArgumentRules();
    ArgumentRules* includeTaxaArgRules2         = new ArgumentRules();
    ArgumentRules* isSequenceMissingArgRules    = new ArgumentRules();
    ArgumentRules* ntaxaArgRules                = new ArgumentRules();
    ArgumentRules* percentageMissingArgRules    = new ArgumentRules();
    ArgumentRules* showdataArgRules             = new ArgumentRules();
    ArgumentRules* removeTaxaArgRules           = new ArgumentRules();
    ArgumentRules* removeTaxaArgRules2          = new ArgumentRules();
    ArgumentRules* setHomeologPhaseArgRules     = new ArgumentRules();
    ArgumentRules* setTaxonNameArgRules         = new ArgumentRules();
    ArgumentRules* setTaxonObjectArgRules       = new ArgumentRules();
    ArgumentRules* taxaArgRules                 = new ArgumentRules();
    ArgumentRules* taxonIndexArgRules           = new ArgumentRules();
    
    std::vector<TypeSpec> taxon_types;
    taxon_types.push_back( RlString::getClassTypeSpec() );
    taxon_types.push_back( Taxon::getClassTypeSpec() );
    taxon_types.push_back( ModelVector<RlString>::getClassTypeSpec() );
    taxon_types.push_back( ModelVector<Taxon>::getClassTypeSpec() );
    addTaxonArgRules->push_back(            new ArgumentRule("taxon" , taxon_types, "The name(s) of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    std::vector<TypeSpec> exclude_taxa_types;
    exclude_taxa_types.push_back( RlString::getClassTypeSpec() );
    exclude_taxa_types.push_back( Taxon::getClassTypeSpec() );
    excludeTaxaArgRules->push_back(         new ArgumentRule("taxon" , exclude_taxa_types, "The name(s) of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    std::vector<TypeSpec> exclude_taxa_types_2;
    exclude_taxa_types_2.push_back( ModelVector<RlString>::getClassTypeSpec() );
    exclude_taxa_types_2.push_back( ModelVector<Taxon>::getClassTypeSpec() );
    excludeTaxaArgRules2->push_back(        new ArgumentRule("taxa" , exclude_taxa_types_2, "The name(s) of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    includeTaxaArgRules->push_back(         new ArgumentRule("name" , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    includeTaxaArgRules2->push_back(        new ArgumentRule("names" , ModelVector<RlString>::getClassTypeSpec(), "The names of the taxa.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    isSequenceMissingArgRules->push_back(   new ArgumentRule("name" , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    percentageMissingArgRules->push_back(   new ArgumentRule("name" , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    removeTaxaArgRules->push_back(          new ArgumentRule("name" , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    removeTaxaArgRules2->push_back(         new ArgumentRule("names" , ModelVector<RlString>::getClassTypeSpec(), "The names of the taxa.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setHomeologPhaseArgRules->push_back(    new ArgumentRule("data_name"  , RlString::getClassTypeSpec(), "The name used in the character alignment.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setHomeologPhaseArgRules->push_back(    new ArgumentRule("tip_name"   , RlString::getClassTypeSpec(), "The tip name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setTaxonNameArgRules->push_back(        new ArgumentRule("current"    , RlString::getClassTypeSpec(), "The old name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setTaxonNameArgRules->push_back(        new ArgumentRule("new"        , RlString::getClassTypeSpec(), "The new name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setTaxonObjectArgRules->push_back(      new ArgumentRule("current"    , RlString::getClassTypeSpec(), "The old name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setTaxonObjectArgRules->push_back(      new ArgumentRule("new"        , Taxon::getClassTypeSpec(), "The new name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    taxonIndexArgRules->push_back(          new ArgumentRule("name"       , RlString::getClassTypeSpec(), "he name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

    methods.addFunction( new MemberProcedure( "addMissingTaxa",  RlUtils::Void, addTaxonArgRules ) );
    methods.addFunction( new MemberProcedure( "excludeTaxa",  RlUtils::Void, excludeTaxaArgRules ) );
    methods.addFunction( new MemberProcedure( "excludeTaxa", RlUtils::Void, excludeTaxaArgRules2 ) );
    methods.addFunction( new MemberProcedure( "filename", RlString::getClassTypeSpec(), filenameArgRules ) );
    methods.addFunction( new MemberProcedure( "includeTaxa", RlUtils::Void, includeTaxaArgRules ) );
    methods.addFunction( new MemberProcedure( "includeTaxa", RlUtils::Void, includeTaxaArgRules2 ) );
    methods.addFunction( new MemberProcedure( "isSequenceMissing", RlBoolean::getClassTypeSpec(), isSequenceMissingArgRules ) );
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules ) );
    methods.addFunction( new MemberProcedure( "ntaxa", Natural::getClassTypeSpec(), ntaxaArgRules ) );
    methods.addFunction( new MemberProcedure( "percentageMissing", Probability::getClassTypeSpec(), percentageMissingArgRules ) );
    methods.addFunction( new MemberProcedure( "show", RlUtils::Void, showdataArgRules ) );
    methods.addFunction( new MemberProcedure( "removeTaxa", RlUtils::Void, removeTaxaArgRules ) );
    methods.addFunction( new MemberProcedure( "removeTaxa", RlUtils::Void, removeTaxaArgRules2 ) );
    methods.addFunction( new MemberProcedure( "setHomeologPhase", RlUtils::Void, setHomeologPhaseArgRules ) );
    methods.addFunction( new MemberProcedure( "setTaxonName", RlUtils::Void, setTaxonNameArgRules ) );
    methods.addFunction( new MemberProcedure( "setTaxonObject", RlUtils::Void, setTaxonObjectArgRules ) );
    methods.addFunction( new MemberProcedure( "taxa", ModelVector<Taxon>::getClassTypeSpec(), taxaArgRules ) );
    methods.addFunction( new MemberProcedure( "taxonIndex", Natural::getClassTypeSpec(), taxonIndexArgRules ) );
    
    // Add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );
    
    return methods;
}

/* Map calls to member methods */
RevPtr<RevVariable> AbstractCharacterData::executeCharacterDataMethod(std::string const &name, const std::vector<Argument> &args, bool &found,  RevBayesCore::AbstractCharacterData *charDataObject)
{
    
    if (name == "addMissingTaxa")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( RlString::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const RlString&>( argument ).getValue();
            charDataObject->addMissingTaxon( n );
        }
        else if ( argument.isType( ModelVector<RlString>::getClassTypeSpec() ) )
        {
            const ModelVector<RlString>& n = static_cast<const ModelVector<RlString>&>( argument );
            for ( size_t i=0; i<n.size(); ++i )
            {
                charDataObject->addMissingTaxon( n[i] );
            }
        }
        else if ( argument.isType( Taxon::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const Taxon&>( argument ).getValue().getName();
            charDataObject->addMissingTaxon( n );
        }
        else if ( argument.isType( ModelVector<Taxon>::getClassTypeSpec() ) )
        {
            const ModelVector<Taxon>& n = static_cast<const ModelVector<Taxon>&>( argument );
            for ( size_t i=0; i<n.size(); ++i )
            {
                charDataObject->addMissingTaxon( n[i].getName() );
            }
        }
        
        return NULL;
    }
    else if (name == "chartype")
    {
        found = true;
        
        return new RevVariable( new RlString( charDataObject->getDataType() ) );
    }
    else if ( name == "removeTaxa" )
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( RlString::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const RlString&>( argument ).getValue();
            charDataObject->deleteTaxon( n );
        }
        else if ( argument.isType( ModelVector<RlString>::getClassTypeSpec() ) )
        {
            const ModelVector<RlString>& x = static_cast<const ModelVector<RlString>&>( argument );
            RevBayesCore::AbstractCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                v.deleteTaxon( x[i] );
            }
        }
        return NULL;
    }
    else if (name == "excludeTaxa" )
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( RlString::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const RlString&>( argument ).getValue();
            charDataObject->excludeTaxon( n );
        }
        else if ( argument.isType( Taxon::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const Taxon&>( argument ).getValue().getSpeciesName();
            charDataObject->excludeTaxon( n );
        }
        else if ( argument.isType( ModelVector<RlString>::getClassTypeSpec() ) )
        {
            const ModelVector<RlString>& x = static_cast<const ModelVector<RlString>&>( argument );
            RevBayesCore::AbstractCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                v.excludeTaxon( x[i] );
            }
        }
        else if ( argument.isType( ModelVector<Taxon>::getClassTypeSpec() ) )
        {
            const ModelVector<Taxon>& x = static_cast<const ModelVector<Taxon>&>( argument );
            RevBayesCore::AbstractCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                v.excludeTaxon( x[i].getSpeciesName() );
            }
        }
        return NULL;
    }
    else if (name == "includeTaxa")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( RlString::getClassTypeSpec() ) )
        {
            const std::string &n = static_cast<const RlString&>( argument ).getValue();
            // remember that we internally store the character indeces from 0 to n-1
            // but externally represent it as 1 to n
            charDataObject->includeTaxon( n );
        }
        else if ( argument.isType( ModelVector<RlString>::getClassTypeSpec() ) )
        {
            const ModelVector<RlString>& x = static_cast<const ModelVector<RlString>&>( argument );
            RevBayesCore::AbstractCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                v.includeTaxon( x[i] );
            }
        }
        return NULL;
    }
    else if (name == "filename")
    {
        found = true;
        
        const std::string n = charDataObject->getFilename().string();
        
        return new RevVariable( new RlString(n) );
    }
    else if (name == "isSequenceMissing")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        const std::string &n = static_cast<const RlString&>( argument ).getValue();
        
        bool tf = charDataObject->isSequenceMissing( n );

        return new RevVariable( new RlBoolean(tf) );
    }
    else if (name == "names" || name == "taxa")
    {
        found = true;
        
        ModelVector<Taxon> *n = new ModelVector<Taxon>();
        for (size_t i = 0; i < charDataObject->getNumberOfTaxa(); ++i)
        {
            n->push_back( charDataObject->getTaxon( i ) );
        }
        
        return new RevVariable( n );
    }
    else if (name == "ntaxa") 
    {
        found = true;
        
        int n = (int)charDataObject->getNumberOfTaxa();
        
        return new RevVariable( new Natural(n) );
    }
    else if (name == "percentageMissing")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        const std::string &n = static_cast<const RlString&>( argument ).getValue();
        
        double p = charDataObject->getPercentageMissing( n );
        
        return new RevVariable( new Probability(p) );
    }
    else if (name == "size") 
    {
        found = true;
        
        int n = (int)charDataObject->getNumberOfTaxa();
        
        return new RevVariable( new Natural(n) );
    }
    else if (name == "setHomeologPhase")
    {
        found = true;
        
        const RevObject& data_name = args[0].getVariable()->getRevObject();
        if ( data_name.isType( RlString::getClassTypeSpec() ) )
        {
            std::string n = std::string( static_cast<const RlString&>( data_name ).getValue() );
            const RevObject& tip_name = args[1].getVariable()->getRevObject();
            if ( tip_name.isType( RlString::getClassTypeSpec() ) )
            {
                std::string name = std::string( static_cast<const RlString&>( tip_name ).getValue() );
                charDataObject->setHomeologPhase( n ,name );
            }
        }
        return NULL;
    }
    else if (name == "setTaxonName")
    {
        found = true;
        
        const RevObject& current = args[0].getVariable()->getRevObject();
        if ( current.isType( RlString::getClassTypeSpec() ) )
        {
            std::string n = std::string( static_cast<const RlString&>( current ).getValue() );
            const RevObject& new_name = args[1].getVariable()->getRevObject();
            if ( new_name.isType( RlString::getClassTypeSpec() ) )
            {
                std::string name = std::string( static_cast<const RlString&>( new_name ).getValue() );
                charDataObject->setTaxonName( n ,name );
               // std::cout << "new name: "<< dagNode->getValue().getTaxonData( n ).getTaxonName() << std::endl;
            }
        }
        return NULL;
    }
    else if (name == "setTaxonObject")
    {
        found = true;
        
        const RevObject& current = args[0].getVariable()->getRevObject();
        if ( current.isType( RlString::getClassTypeSpec() ) )
        {
            std::string n = std::string( static_cast<const RlString&>( current ).getValue() );
            const RevObject& new_taxon_object = args[1].getVariable()->getRevObject();
            if ( new_taxon_object.isType( Taxon::getClassTypeSpec() ) )
            {
                const RevBayesCore::Taxon &t = static_cast<const Taxon&>( new_taxon_object ).getValue();
                charDataObject->setTaxonObject( n ,t );
                // std::cout << "new name: "<< dagNode->getValue().getTaxonData( n ).getTaxonName() << std::endl;
            }
        }
        return NULL;
    }
    else if (name == "show")
    {
        found = true;
        
        charDataObject->show(std::cout);
        
        return NULL;
    }
    else if (name == "taxonIndex")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        const std::string &n = static_cast<const RlString&>( argument ).getValue();
        int index = (int)charDataObject->getIndexOfTaxon( n ) + 1;
        
        return new RevVariable( new Natural(index) );
    }
    
    
    // not found a matching method
    found = false;
    return NULL;
}


