#ifndef CharacterTranslator_H
#define CharacterTranslator_H

#include <iosfwd>

#include "AbstractDiscreteTaxonData.h"
#include "DiscreteTaxonData.h"

namespace RevBayesCore {
class AminoAcidState;
class CodonState;
class DnaState;
class RnaState;
class DoubletState;
    
    namespace CharacterTranslator {
    
    
        AbstractDiscreteTaxonData*              translateCharacters(const AbstractDiscreteTaxonData &d, const std::string& type);             
    
        DiscreteTaxonData<CodonState>*          translateToCodon(const AbstractDiscreteTaxonData &d);
        DiscreteTaxonData<CodonState>*          translateToCodonFromDna(const DiscreteTaxonData<DnaState> &d);
        DiscreteTaxonData<CodonState>*          translateToCodonFromRna(const DiscreteTaxonData<RnaState> &d);

        DiscreteTaxonData<CodonState>*          translateToDoublet(const AbstractDiscreteTaxonData &d);
        DiscreteTaxonData<DoubletState>*        translateToDoubletFromDna(const DiscreteTaxonData<DnaState> &d);
        DiscreteTaxonData<DoubletState>*        translateToDoubletFromRna(const DiscreteTaxonData<RnaState> &d);

        DiscreteTaxonData<DnaState>*            translateToDna(const AbstractDiscreteTaxonData &d);
//        DiscreteTaxonData<DnaState>*            translateToDna(const DiscreteTaxonData<AminoAcidState> &d);                                   
//        DiscreteTaxonData<DnaState>*            translateToDna(const DiscreteTaxonData<CodonState> &d);                                       
//        DiscreteTaxonData<DnaState>*            translateToDna(const DiscreteTaxonData<DnaState> &d);                                         
//        DiscreteTaxonData<DnaState>*            translateToDna(const DiscreteTaxonData<RnaState> &d);                                         
        
        DiscreteTaxonData<RnaState>*            translateToRna(const AbstractDiscreteTaxonData &d);                                           
//        DiscreteTaxonData<RnaState>*            translateToRna(const DiscreteTaxonData<AminoAcidState> &d);                                   

        
        DiscreteTaxonData<AminoAcidState>*      translateToAa(const AbstractDiscreteTaxonData &d);                                            
//        DiscreteTaxonData<AminoAcidState>*      translateToAa(const DiscreteTaxonData<AminoAcidState> &d);                                    
        DiscreteTaxonData<AminoAcidState>*      translateToAaFromDna(const DiscreteTaxonData<DnaState> &d);

        
        
    }
    
}

#endif
