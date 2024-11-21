#include "CharacterTranslator.h"

#include <cstddef>
#include <string>

#include "RbException.h"
#include "StringUtilities.h"
#include "AminoAcidState.h"
#include "Cloneable.h"
#include "CodonState.h"
#include "DoubletState.h"
#include "DiscreteCharacterState.h"
#include "DnaState.h"
#include "NaturalNumbersState.h"
#include "RnaState.h"



using namespace RevBayesCore;

RevBayesCore::AbstractDiscreteTaxonData* CharacterTranslator::translateCharacters(const RevBayesCore::AbstractDiscreteTaxonData &d, const std::string &type)
{
    AbstractDiscreteTaxonData *trans_taxon_data = NULL;

//    options.push_back( "DNA" );
//    options.push_back( "RNA" );
//    options.push_back( "AA" );
//    options.push_back( "PoMo" );
//    options.push_back( "Protein" );
//    options.push_back( "Standard" );
//    options.push_back( "NaturalNumbers" );
//    options.push_back( "Binary" );
    
    if ( type == "DNA" )
    {
        trans_taxon_data = translateToDna( d );
    }
    else if ( type == "RNA" )
    {
        trans_taxon_data = translateToRna( d );
    }
    else if ( type == "AA" )
    {
        if ( d.getCharacter(0).getDataType() == "DNA" )
        {
            trans_taxon_data = translateToAaFromDna( dynamic_cast< const DiscreteTaxonData<DnaState>& >(d) );
        }
        else
        {
            trans_taxon_data = translateToAa( d );
        }
    }
    else if ( type == "Protein" )
    {
        if ( d.getCharacter(0).getDataType() == "DNA" )
        {
            trans_taxon_data = translateToAaFromDna( dynamic_cast< const DiscreteTaxonData<DnaState>& >(d) );
        }
        else
        {
            trans_taxon_data = translateToAa( d );
        }
    }
    else if ( type == "Codon" )
    {
        if ( d.getCharacter(0).getDataType() == "DNA" )
        {
            trans_taxon_data = translateToCodonFromDna( dynamic_cast< const DiscreteTaxonData<DnaState>& >(d) );
        }
        else if ( d.getCharacter(0).getDataType() == "RNA" )
        {
            trans_taxon_data = translateToCodonFromRna( dynamic_cast< const DiscreteTaxonData<RnaState>& >(d) );
        }
        else
        {
            trans_taxon_data = translateToCodon( d );
        }
    }
    else if ( type == "Doublet")
    {
        if ( d.getCharacter(0).getDataType() == "DNA" )
        {
            trans_taxon_data = translateToDoubletFromDna( dynamic_cast< const DiscreteTaxonData<DnaState>& >(d) );
        }
        else if ( d.getCharacter(0).getDataType() == "RNA" )
        {
            trans_taxon_data = translateToDoubletFromRna( dynamic_cast< const DiscreteTaxonData<RnaState>& >(d) );
        }
        else
        {
            trans_taxon_data = translateToDoublet( d );
        }
    }
    else
    {
        throw RbException() << "Cannot translate character data object of type \""<< d.getCharacter(0).getDataType() <<"\" into type \"" << type << "\"" ;
    }
    
    return trans_taxon_data;
}


DiscreteTaxonData<DnaState>* CharacterTranslator::translateToDna(const AbstractDiscreteTaxonData &d)
{
    throw RbException("Cannot translate character data object of type \""+ d.getCharacter(0).getDataType() +"\" into type \"DNA\"" );
    
    return NULL;
}


DiscreteTaxonData<RnaState>* CharacterTranslator::translateToRna(const AbstractDiscreteTaxonData &d)
{
    throw RbException("Cannot translate character data object of type \""+ d.getCharacter(0).getDataType() +"\" into type \"RNA\"" );
    
    return NULL;
}


DiscreteTaxonData<AminoAcidState>* CharacterTranslator::translateToAa(const AbstractDiscreteTaxonData &d)
{
    throw RbException("Cannot translate character data object of type \""+ d.getCharacter(0).getDataType() +"\" into type \"AA\"" );
    
    return NULL;
}


DiscreteTaxonData<CodonState>* CharacterTranslator::translateToCodon(const AbstractDiscreteTaxonData &d)
{
    throw RbException("Cannot translate character data object of type \""+ d.getCharacter(0).getDataType() +"\" into type \"Codon\"" );
    
    return NULL;
}


DiscreteTaxonData<CodonState>* CharacterTranslator::translateToDoublet(const AbstractDiscreteTaxonData &d)
{
    throw RbException("Cannot translate character data object of type \""+ d.getCharacter(0).getDataType() +"\" into type \"Doublet\"" );

    return NULL;
}


DiscreteTaxonData<DoubletState>* CharacterTranslator::translateToDoubletFromDna(const DiscreteTaxonData<DnaState> &d)
{
    size_t length = d.getNumberOfCharacters();

    DiscreteTaxonData<DoubletState> *trans_taxon_data = new DiscreteTaxonData<DoubletState>( d.getTaxon() );

    for (size_t i=0; i<(length-1); i+=2)
    {
        std::string doublet_string = d.getCharacter(i).getStringValue();
        doublet_string += d.getCharacter(i+1).getStringValue();
        DoubletState cs = DoubletState( doublet_string );

        trans_taxon_data->addCharacter( cs );
    }

    return trans_taxon_data;
}


DiscreteTaxonData<DoubletState>* CharacterTranslator::translateToDoubletFromRna(const DiscreteTaxonData<RnaState> &d)
{
    size_t length = d.getNumberOfCharacters();

    DiscreteTaxonData<DoubletState> *trans_taxon_data = new DiscreteTaxonData<DoubletState>( d.getTaxon() );

    for (size_t i=0; i<(length-1); i+=2)
    {
        std::string doublet_string = d.getCharacter(i).getStringValue();
        doublet_string += d.getCharacter(i+1).getStringValue();
        StringUtilities::replaceSubstring(doublet_string, "U", "T");
        DoubletState cs = DoubletState( doublet_string );

        trans_taxon_data->addCharacter( cs );
    }

    return trans_taxon_data;
}


DiscreteTaxonData<CodonState>* CharacterTranslator::translateToCodonFromDna(const DiscreteTaxonData<DnaState> &d)
{
    size_t length = d.getNumberOfCharacters();
    
    DiscreteTaxonData<CodonState> *trans_taxon_data = new DiscreteTaxonData<CodonState>( d.getTaxon() );
    
    for (size_t i=0; i<(length-2); i+=3)
    {
        std::string codon_string = d.getCharacter(i).getStringValue();
        codon_string += d.getCharacter(i+1).getStringValue();
        codon_string += d.getCharacter(i+2).getStringValue();
        CodonState cs = CodonState( codon_string );
        
        trans_taxon_data->addCharacter( cs );
    }
    
    return trans_taxon_data;
}


DiscreteTaxonData<CodonState>* CharacterTranslator::translateToCodonFromRna(const DiscreteTaxonData<RnaState> &d)
{
    size_t length = d.getNumberOfCharacters();
    
    DiscreteTaxonData<CodonState> *trans_taxon_data = new DiscreteTaxonData<CodonState>( d.getTaxon() );
    
    for (size_t i=0; i<(length-2); i+=3)
    {
        std::string codon_string = d.getCharacter(i).getStringValue();
        codon_string += d.getCharacter(i+1).getStringValue();
        codon_string += d.getCharacter(i+2).getStringValue();
        StringUtilities::replaceSubstring(codon_string, "U", "T");
        CodonState cs = CodonState( codon_string );
        
        trans_taxon_data->addCharacter( cs );
    }
    
    return trans_taxon_data;
}


DiscreteTaxonData<AminoAcidState>* CharacterTranslator::translateToAaFromDna(const DiscreteTaxonData<DnaState> &d)
{
    size_t length = d.getNumberOfCharacters();
    
    DiscreteTaxonData<AminoAcidState> *trans_taxon_data = new DiscreteTaxonData<AminoAcidState>( d.getTaxon() );
    
    for (size_t i=0; i<(length-2); i+=3)
    {
        std::string codon_string = d.getCharacter(i).getStringValue();
        codon_string += d.getCharacter(i+1).getStringValue();
        codon_string += d.getCharacter(i+2).getStringValue();
        CodonState cs = CodonState( codon_string );
        
        trans_taxon_data->addCharacter( cs.getAminoAcidState() );
    }
    
    return trans_taxon_data;
}


