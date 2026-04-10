#ifndef CharTypeApply_H
#define CharTypeApply_H

#include "BinaryState.h"
#include "NaturalNumbersState.h"
#include "DnaState.h"
#include "RnaState.h"
#include "CodonState.h"
#include "DoubletState.h"
#include "StandardState.h"
#include "AminoAcidState.h"
#include "PoMoState.h"

template <typename F>
void apply_to_character_type(const F& f, const std::string& character_type)
{
    namespace Core = RevBayesCore;

    if (character_type == "AA" || character_type == "Protein")
    {
        // for complicated reasons we can't just write f<Core::AminoAcidState>()
        f.template operator()<Core::AminoAcidState>();
    }
    else if (character_type == "DNA")
    {
        f.template operator()<Core::DnaState>();
    }
    else if (character_type == "Codon")
    {
        f.template operator()<Core::CodonState>();
    }
    else if (character_type == "Doublet")
    {
        f.template operator()<Core::DoubletState>();
    }
    else if (character_type == "NaturalNumbers")
    {
        f.template operator()<Core::NaturalNumbersState>();
    }
    else if (character_type == "PoMo")
    {
        f.template operator()<Core::PoMoState>();
    }
    else if (character_type == "RNA")
    {
        f.template operator()<Core::RnaState>();
    }
    else if (character_type == "Standard")
    {
        f.template operator()<Core::StandardState>();
    }
    else if (character_type == "Binary" || character_type == "Restriction")
    {
        f.template operator()<Core::BinaryState>();
    }
    else
    {
        throw RbException()<<"Incorrect character type '"<<character_type<<"' specified. Valid options are: DNA, RNA, AA or Protein, Codon, Doublet, NaturalNumbers, PoMo, Standard, Binary or Restriction";
    }
}

#endif
