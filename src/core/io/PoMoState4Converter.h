#ifndef PoMoState4Converter_H
#define PoMoState4Converter_H

#include <map>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousDiscreteCharacterData.h"
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
    class PoMoState4Converter {
        
    public:
        PoMoState4Converter();
        
        HomologousDiscreteCharacterData<PoMoState4>*  convertData(const AbstractHomologousDiscreteCharacterData &d, const unsigned int virtualPopulationSize, const std::map<std::string, std::string> sequenceNameToSpeciesName);
        
    private:
        PoMoState4* convertCounts(std::vector<double> &counts,
                                        const unsigned int virtualPopulationSize,
                                        std::vector< std::vector<double> > &frequencies);
        
    };
    
}


#endif /* defined(PoMoState4Converter_H) */
