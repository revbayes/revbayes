#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "DelimitedCharacterDataReader.h"
#include "HomologousDiscreteCharacterData.h"
#include "Ellipsis.h"
#include "Func_readCharacterDataDelimited.h"
#include "NaturalNumbersState.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlContinuousCharacterData.h"
#include "RlNaturalNumbersState.h"
#include "RlStandardState.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readCharacterDataDelimited* Func_readCharacterDataDelimited::clone( void ) const
{
    
    return new Func_readCharacterDataDelimited( *this );
}


std::string Func_readCharacterDataDelimited::bitToState(const std::string &s)
{
    
    std::stringstream ss;
    char* ptr;
    long parsed = strtol(s.c_str(), &ptr, 2);
    
    if (parsed > RbConstants::Integer::max)
    {
        throw RbException("ERROR: readTSVBitsetData token " + s + " too large to store as NaturalNumber");
    }
    
    if (s.find_first_not_of("01") != std::string::npos)
    {
        throw RbException("ERROR: readTSVBitsetData token " + s + " contains non-binary characters");
    }
    
    ss << parsed;
    
    return ss.str();
}


/** Execute function */
RevPtr<RevVariable> Func_readCharacterDataDelimited::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string&  fn      = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  dt      = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string&  lab     = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();
    const std::string&  del     = static_cast<const RlString&>( args[3].getVariable()->getRevObject() ).getValue();
    bool                header  = static_cast<const RlBoolean&>( args[4].getVariable()->getRevObject() ).getValue();
    size_t lines_to_skip = ( header == true ? 1 : 0 );
   
    if (lab.compare(std::string("")) == 0)
    {
        throw RbException("You must enter a value for the state labels.");
    }

    if (dt == "NaturalNumbers")
    {
        
        // setup a matrix of NaturalNumbers states
        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::NaturalNumbersState> *coreStates = new RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::NaturalNumbersState>();
        
        // get data from file
        RevBayesCore::DelimitedCharacterDataReader* tsv_data = new RevBayesCore::DelimitedCharacterDataReader(fn, del[0], lines_to_skip);

        int max = StringUtilities::asIntegerNumber( lab );
        
        // loop through data and get each NaturalNumbers value
        for (size_t i = 0; i < tsv_data->getData().size(); ++i)
        {
            
            // now put core state into DiscreteTaxonData
            RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState> coreSeq = RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState>(tsv_data->getNames()[i]);

            // get count from data
            const std::vector<std::string> &data = tsv_data->getData()[i];
            
            for (size_t j= 0; j < data.size(); ++j)
            {
                // make the core state
                RevBayesCore::NaturalNumbersState coreState = RevBayesCore::NaturalNumbersState( data[j], max );
                
                coreSeq.addCharacter( coreState );
            }
            
            // add DiscreteTaxonData to the matrix of NaturalNumbers states
            coreStates->addTaxonData( coreSeq );
        
        }
        
        // put coreNaturalNumbers matrix into rev language level matrix
//        HomologousDiscreteCharacterData<RevLanguage::NaturalNumbersState> *rlCoreStates = new HomologousDiscreteCharacterData<RevLanguage::NaturalNumbersState>( coreStates );
        AbstractHomologousDiscreteCharacterData *rlCoreStates = new AbstractHomologousDiscreteCharacterData( coreStates );
        
        return new RevVariable( rlCoreStates );
        
    }
    else if (dt == "Bitset")
    {
    
        // setup a matrix of NaturalNumbers states
        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::NaturalNumbersState> *coreStates = new RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::NaturalNumbersState>();
        
        // get data from file
        RevBayesCore::DelimitedCharacterDataReader* tsv_data = new RevBayesCore::DelimitedCharacterDataReader(fn, del[0], lines_to_skip);
        
        // loop through data and get each NaturalNumbers value
        for (size_t i = 0; i < tsv_data->getData().size(); ++i)
        {
            
            // now put core state into DiscreteTaxonData
            RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState> coreSeq = RevBayesCore::DiscreteTaxonData<RevBayesCore::NaturalNumbersState>(tsv_data->getNames()[i]);
            
            // get count from data
            const std::vector<std::string> &data = tsv_data->getData()[i];
            
            for (size_t j= 0; j < data.size(); ++j)
            {
                // make the core state
                std::string d = data[j];
                std::reverse(d.begin(), d.end());
                d = bitToState( d );
                RevBayesCore::NaturalNumbersState coreState = RevBayesCore::NaturalNumbersState( d, 1 );
                
                coreSeq.addCharacter( coreState );
            }
            
            // add DiscreteTaxonData to the matrix of NaturalNumbers states
            coreStates->addTaxonData( coreSeq );
            
        }
        
        // put coreNaturalNumbers matrix into rev language level matrix
//        HomologousDiscreteCharacterData<RevLanguage::NaturalNumbersState> *rlCoreStates = new HomologousDiscreteCharacterData<RevLanguage::NaturalNumbersState>( coreStates );
        AbstractHomologousDiscreteCharacterData *rlCoreStates = new AbstractHomologousDiscreteCharacterData( coreStates );
        
        return new RevVariable( rlCoreStates );
        
    }
    else if (dt == "Standard")
    {
        
        // setup a matrix of NaturalNumbers states
        RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState> *coreStates = new RevBayesCore::HomologousDiscreteCharacterData<RevBayesCore::StandardState>();
        
        // get data from file
        RevBayesCore::DelimitedCharacterDataReader* tsv_data = new RevBayesCore::DelimitedCharacterDataReader(fn, del[0], lines_to_skip);
        
        // loop through data and get each NaturalNumbers value
        for (size_t i = 0; i < tsv_data->getData().size(); ++i)
        {
            
            // now put core state into DiscreteTaxonData
            RevBayesCore::DiscreteTaxonData<RevBayesCore::StandardState> coreSeq = RevBayesCore::DiscreteTaxonData<RevBayesCore::StandardState>(tsv_data->getNames()[i]);
            
            // get count from data
            const std::vector<std::string> &data = tsv_data->getData()[i];
            
            for (size_t j= 0; j < data.size(); ++j)
            {
                // make the core state
                RevBayesCore::StandardState coreState = RevBayesCore::StandardState( data[j], lab );
                
                coreSeq.addCharacter( coreState );
            }
            
            // add DiscreteTaxonData to the matrix of NaturalNumbers states
            coreStates->addTaxonData( coreSeq );
            
        }
        
        // put coreNaturalNumbers matrix into rev language level matrix
//        HomologousDiscreteCharacterData<RevLanguage::StandardState> *rlCoreStates = new HomologousDiscreteCharacterData<RevLanguage::StandardState>( coreStates );
        AbstractHomologousDiscreteCharacterData *rlCoreStates = new AbstractHomologousDiscreteCharacterData( coreStates );
        
        return new RevVariable( rlCoreStates );
        
    }
    else if (dt == "Continuous")
    {
        
        // setup a matrix of NaturalNumbers states
        RevBayesCore::ContinuousCharacterData *coreStates = new RevBayesCore::ContinuousCharacterData();
        
        // get data from file
        RevBayesCore::DelimitedCharacterDataReader* tsv_data = new RevBayesCore::DelimitedCharacterDataReader(fn, del[0], lines_to_skip);
        
        // loop through data and get each NaturalNumbers value
        for (size_t i = 0; i < tsv_data->getData().size(); ++i)
        {
            
            // now put core state into DiscreteTaxonData
            RevBayesCore::ContinuousTaxonData coreSeq = RevBayesCore::ContinuousTaxonData(tsv_data->getNames()[i]);
            
            // get count from data
            const std::vector<std::string> &data = tsv_data->getData()[i];
            
            for (size_t j= 0; j < data.size(); ++j)
            {
                // make the core state
                double coreState = atof( data[j].c_str() );
                
                coreSeq.addCharacter( coreState );
            }
            
            // add DiscreteTaxonData to the matrix of NaturalNumbers states
            coreStates->addTaxonData( coreSeq );
            
        }
        
        // put coreNaturalNumbers matrix into rev language level matrix
        ContinuousCharacterData *rlCoreStates = new ContinuousCharacterData( coreStates );
        
        return new RevVariable( rlCoreStates );
        
    }
    else
    {
        
        throw RbException( "Invalid data type. Valid data types are: NaturalNumbers|Bitset|Standard|Continuous" );
        
    }
}


/** Get argument rules */
const ArgumentRules& Func_readCharacterDataDelimited::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file",      RlString::getClassTypeSpec(), "The name of the file to read in.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        std::vector<std::string> type_options;
        type_options.push_back( "NaturalNumbers" );
        type_options.push_back( "Bitset" );
        type_options.push_back( "Standard" );
        type_options.push_back( "Continuous" );
        argumentRules.push_back( new OptionRule( "type", new RlString("NaturalNumbers"), type_options, "The type of data." ) );
        argumentRules.push_back( new ArgumentRule( "stateLabels", RlString::getClassTypeSpec(), "The state labels (for standard states) or max number for NaturalNumbers.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "" ) ) );
        argumentRules.push_back( new ArgumentRule( "delimiter", RlString::getClassTypeSpec(), "The delimiter between columns.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "\t" ) ) );
        argumentRules.push_back( new ArgumentRule( "headers", RlBoolean::getClassTypeSpec(), "Has this file a header line?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        rules_set = true;
        
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readCharacterDataDelimited::getClassType(void)
{
    
    static std::string rev_type = "Func_readCharacterDataDelimited";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readCharacterDataDelimited::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readCharacterDataDelimited::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readCharacterDataDelimited";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readCharacterDataDelimited::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readCharacterDataDelimited::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = NaturalNumbersState::getClassTypeSpec();
    return return_typeSpec;
}




