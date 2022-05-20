#ifndef PoMoStateConverter_H
#define PoMoStateConverter_H

#include <map>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousDiscreteCharacterData.h"
#include "PoMoState.h"


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
        
        HomologousDiscreteCharacterData<PoMoState>*     convertData2(const AbstractHomologousDiscreteCharacterData &d, size_t virtual_ps, const RbVector<Taxon>& taxa);
        HomologousDiscreteCharacterData<PoMoState>*     convertData4(const AbstractHomologousDiscreteCharacterData &d, size_t virtual_ps, const RbVector<Taxon>& taxa);
        
    private:
        void                                            createSpeciesToSampleNamesMap( std::map<std::string, std::vector<std::string> >& s, const RbVector<Taxon>& taxa );

    };
    
}


#endif /* defined(PoMoStateConverter_H) */
