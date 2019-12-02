#include "RlNonHomologousCharacterData.h"
#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RlString.h"


using namespace RevLanguage;

NonHomologousCharacterData::NonHomologousCharacterData( void ) : AbstractCharacterData( )
{
    
}


MethodTable NonHomologousCharacterData::getCharacterDataMethods( void ) const
{
    
    MethodTable methods = AbstractCharacterData::getCharacterDataMethods();
    
    ArgumentRules* ncharArgRules                = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure("nchar", Natural::getClassTypeSpec(),   ncharArgRules        ) );
    
    
    return methods;
}


/* Map calls to member methods */
RevPtr<RevVariable> NonHomologousCharacterData::executeCharacterDataMethod(std::string const &name, const std::vector<Argument> &args, bool &found, RevBayesCore::AbstractCharacterData *charDataObject)
{
    RevPtr<RevVariable> retVal = AbstractCharacterData::executeCharacterDataMethod(name, args, found, charDataObject);
 
//    RevBayesCore::NonHomologousCharacterData *charDataObject = &getValue();

//    if (name == "nchar")
//    {
//        found = true;
//        
//        Natural *numChar = new Natural( charDataObject->getNumberOfIncludedCharacters() );
//        
//        return new RevVariable( numChar );
//    }
    
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


