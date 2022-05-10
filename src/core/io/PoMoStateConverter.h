#ifndef PoMoStateConverter_H
#define PoMoStateConverter_H

#include <map>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousDiscreteCharacterData.h"
#include "PoMoState.h"
#include "PoMoState4.h"


namespace RevBayesCore {
    
    /**
     * This class enables conversion of a DNA matrix into a PoMoState4 matrix.
     * The PoMoState4 matrix can then be analyzed with the PoMo model.
     *
     * This class currently has only one functionality,
     * to convert the DNA data into PoMoState4 data.
     *
     */
    class PoMoStateConverter {
        
    public:
        PoMoStateConverter();
        
        HomologousDiscreteCharacterData<PoMoState>*   convertData2(const AbstractHomologousDiscreteCharacterData &d, size_t virtual_ps, const RbVector<Taxon>& taxa);
        HomologousDiscreteCharacterData<PoMoState4>*  convertData4(const AbstractHomologousDiscreteCharacterData &d, size_t virtual_ps, const std::map<std::string, std::string> sequenceNameToSpeciesName);
        
    private:
        PoMoState4*                                   convertCounts4(std::vector<double> &counts, size_t virtualPopulationSize, std::vector< std::vector<double> > &frequencies);

    };
    
}


#endif /* defined(PoMoStateConverter_H) */
