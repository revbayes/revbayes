#include "RlHomologousCharacterData.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "HomologousCharacterData.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

HomologousCharacterData::HomologousCharacterData( void ) : AbstractCharacterData(  )
{
    
}


MethodTable HomologousCharacterData::getCharacterDataMethods( void ) const
{
    
    MethodTable methods = AbstractCharacterData::getCharacterDataMethods();
    
    ArgumentRules* ncharArgRules                = new ArgumentRules();
    ArgumentRules* excludeallArgRules           = new ArgumentRules();
    ArgumentRules* includeallArgRules           = new ArgumentRules();
    ArgumentRules* excludecharArgRules          = new ArgumentRules();
    ArgumentRules* excludecharArgRules2         = new ArgumentRules();
    ArgumentRules* includecharArgRules          = new ArgumentRules();
    ArgumentRules* includecharArgRules2         = new ArgumentRules();
    ArgumentRules* isResolvedArgRules           = new ArgumentRules();
    ArgumentRules* getIncludedCharacterIndicesArgRules = new ArgumentRules();

    
    excludecharArgRules->push_back(         new ArgumentRule("pos"        , Natural::getClassTypeSpec()              , "The position of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    excludecharArgRules2->push_back(        new ArgumentRule(""           , ModelVector<Natural>::getClassTypeSpec() , "A vector of character positions.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    includecharArgRules->push_back(         new ArgumentRule(""           , Natural::getClassTypeSpec()              , "The position of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    includecharArgRules2->push_back(        new ArgumentRule(""           , ModelVector<Natural>::getClassTypeSpec() , "A vector of character positions.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    isResolvedArgRules->push_back(          new ArgumentRule("taxonIndex" , Natural::getClassTypeSpec()              , "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    isResolvedArgRules->push_back(          new ArgumentRule("charIndex"  , Natural::getClassTypeSpec()              , "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    
    methods.addFunction( new MemberProcedure( "nchar", Natural::getClassTypeSpec(),         ncharArgRules        ) );
    methods.addFunction( new MemberProcedure( "excludeCharacter", RlUtils::Void,            excludecharArgRules  ) );
    methods.addFunction( new MemberProcedure( "excludeCharacter", RlUtils::Void,            excludecharArgRules2 ) );
    methods.addFunction( new MemberProcedure( "excludeAll", RlUtils::Void,                  excludeallArgRules   ) );
    methods.addFunction( new MemberProcedure( "getIncludedCharacterIndices", ModelVector<Natural>::getClassTypeSpec(), getIncludedCharacterIndicesArgRules   ) );
    methods.addFunction( new MemberProcedure( "includeCharacter", RlUtils::Void,            includecharArgRules  ) );
    methods.addFunction( new MemberProcedure( "includeCharacter", RlUtils::Void,            includecharArgRules2 ) );
    methods.addFunction( new MemberProcedure( "includeAll", RlUtils::Void,                  includeallArgRules   ) );
    methods.addFunction( new MemberProcedure( "isResolved", RlBoolean::getClassTypeSpec(),  isResolvedArgRules   ) );
    
    
    return methods;
}

/* Map calls to member methods */
RevPtr<RevVariable> HomologousCharacterData::executeCharacterDataMethod(std::string const &name, const std::vector<Argument> &args, bool &found, RevBayesCore::HomologousCharacterData *charDataObject)
{
    RevPtr<RevVariable> retVal = AbstractCharacterData::executeCharacterDataMethod(name, args, found, charDataObject);
    
    if (name == "excludeCharacter")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( Natural::getClassTypeSpec() ) )
        {
            size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );
            // remember that we internally store the character indeces from 0 to n-1
            // but externally represent it as 1 to n
            charDataObject->excludeCharacter( n-1 );
        }
        else if ( argument.isType( ModelVector<Natural>::getClassTypeSpec() ) )
        {
            const ModelVector<Natural>& x = static_cast<const ModelVector<Natural>&>( argument );
            RevBayesCore::HomologousCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                // remember that we internally store the character indeces from 0 to n-1
                // but externally represent it as 1 to n
                v.excludeCharacter( size_t(x[i])-1 );
            }
        }
        return NULL;
    }
    else if (name == "excludeAll")
    {
        found = true;
        
        RevBayesCore::HomologousCharacterData &v = *charDataObject;
        size_t nChars = v.getNumberOfCharacters();
        
        for (size_t i = 0; i < nChars; i++)
        {
            v.excludeCharacter(i);
        }
        
        return NULL;
    }
    else if (name == "getIncludedCharacterIndices")
    {
        found = true;
        
        RevBayesCore::HomologousCharacterData &v = *charDataObject;
        size_t nChars = v.getNumberOfCharacters();
        
        RevBayesCore::RbVector<long> indices;
        for (size_t i = 0; i < nChars; i++)
        {
            bool excluded = v.isCharacterExcluded(i);
            if ( excluded == false )
            {
                // remember that we internally store the character indices from 0 to n-1
                // but externally represent it as 1 to n
                indices.push_back( int(i+1) );
            }
        }
        
        return new RevVariable( new ModelVector<Natural>( indices ) );
        
    }
    else if (name == "includeCharacter")
    {
        found = true;
        
        const RevObject& argument = args[0].getVariable()->getRevObject();
        if ( argument.isType( Natural::getClassTypeSpec() ) )
        {
            size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );
            // remember that we internally store the character indices from 0 to n-1
            // but externally represent it as 1 to n
            charDataObject->includeCharacter( n-1 );
        }
        else if ( argument.isType( ModelVector<Natural>::getClassTypeSpec() ) )
        {
            const ModelVector<Natural>& x = static_cast<const ModelVector<Natural>&>( argument );
            RevBayesCore::HomologousCharacterData &v = *charDataObject;
            for ( size_t i=0; i<x.size(); i++ )
            {
                // remember that we internally store the character indices from 0 to n-1
                // but externally represent it as 1 to n
                v.includeCharacter( size_t(x[i])-1 );
            }
        }
        
        return NULL;
    }
    else if (name == "includeAll")
    {
        found = true;
        
        RevBayesCore::HomologousCharacterData &v = *charDataObject;
        size_t nChars = v.getNumberOfCharacters();
        
        for (size_t i = 0; i < nChars; i++)
        {
            v.includeCharacter(i);
        }
        
        return NULL;
        
    }
    else if (name == "isResolved")
    {
        found = true;
        
        const RevObject& argument1 = args[0].getVariable()->getRevObject();
        const RevObject& argument2 = args[1].getVariable()->getRevObject();
        size_t tIdx = size_t( static_cast<const Natural&>( argument1 ).getValue() );
        size_t cIdx = size_t( static_cast<const Natural&>( argument2 ).getValue() );
        // remember that we internally store the character indices from 0 to n-1
        // but externally represent it as 1 to n
        bool tf = charDataObject->isCharacterResolved( tIdx-1, cIdx-1 );
        return new RevVariable( new RlBoolean(tf) );
    }
    else if (name == "nchar")
    {
        found = true;
        
        Natural *numChar = new Natural( charDataObject->getNumberOfIncludedCharacters() );
        
        return new RevVariable( numChar );
    }
    
    //    else if (name == "nexcludedtaxa")
    //    {
    //        int n = (int)deletedTaxa.size();
    //        return RevPtr<RevObject>( new Natural(n) );
    //    }
    //    else if (name == "nexcludedchars")
    //    {
    //        int n = (int)deletedCharacters.size();
    //        return RevPtr<RevObject>( new Natural(n) );
    //    }
    //    else if (name == "nincludedtaxa")
    //    {
    //        int n = (int)(data.getNumberOfTaxa() - deletedTaxa.size());
    //        return RevPtr<RevObject>( new Natural(n) );
    //    }
    //    else if (name == "nincludedchars")
    //    {
    //        int n = (int)(data.getNumberOfCharacters() - deletedCharacters.size());
    //        return RevPtr<RevObject>( new Natural(n) );
    //    }
    //    else if (name == "excludedtaxa")
    //    {
    //        RlVector<std::string, RbString> *excludedTaxa = new RlVector<std::string, RbString>();
    //        for (std::set<size_t>::iterator it = data.deletedTaxa.begin(); it != deletedTaxa.end(); it++)
    //        {
    //            const std::string& tn = getTaxonNameWithIndex(*it);
    //            excludedTaxa->push_back( new RbString( tn ) );
    //        }
    //        return excludedTaxa;
    //    }
    //    else if (name == "excludedchars")
    //    {
    //        RlVector *excludedChars = new RlVector(Natural::getClassTypeSpec());
    //        for (std::set<size_t>::iterator it = deletedCharacters.begin(); it != deletedCharacters.end(); it++)
    //            excludedChars->push_back( new Natural(*it) );
    //        return excludedChars;
    //    }
    //    else if (name == "includedtaxa")
    //    {
    //        RlVector *includedTaxa = new RlVector(RbString::getClassTypeSpec());
    //        for (size_t i=0; i<getNumberOfTaxa(); i++)
    //        {
    //            if ( isTaxonExcluded(i) == false )
    //                includedTaxa->push_back( new RbString( getTaxonNameWithIndex(i) ) );
    //        }
    //        return includedTaxa;
    //    }
    //    else if (name == "includedchars")
    //    {
    //        RlVector *includedChars = new RlVector(Natural::getClassTypeSpec());
    //        for (size_t i=0; i<getNumberOfCharacters(); i++)
    //        {
    //            if ( isCharacterExcluded(i) == false )
    //                includedChars->push_back( new Natural(i+1) );
    //        }
    //        return includedChars;
    //    }
    
    
    // not found a matching method
//    found = false;
    
    return retVal;
}
