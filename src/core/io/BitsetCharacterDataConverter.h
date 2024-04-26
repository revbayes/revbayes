//
//  BitsetCharacterDataConverter.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/1/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef BitsetCharacterDataConverter_h
#define BitsetCharacterDataConverter_h

#include <vector>
#include <map>
#include <cstddef>
#include <iosfwd>

#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "StandardState.h"
#include "DiscreteTaxonData.h"
#include "RbException.h"

namespace RevBayesCore {
    
    class BitsetCharacterDataConverter {
        
    public:
        BitsetCharacterDataConverter(const HomologousDiscreteCharacterData<StandardState>& d, std::string f, size_t ns=0);
//        BitsetCharacterDataConverter(void);
        
        HomologousDiscreteCharacterData<NaturalNumbersState>*  convertData(void);
        // const AbstractHomologousDiscreteCharacterData &d, const size_t int maxRangeSize);
        
    private:
        
        size_t numBitsOn(std::vector<size_t> v);
        void initializeBits(size_t n);
//        void populateRangesRecursively(std::vector<size_t> v, size_t n, size_t idx);
        
        std::map< std::vector<size_t>, size_t >             bitsToStates;
        std::vector<std::vector<size_t> >                   statesToBits;
        std::vector<std::vector<std::vector<size_t> > >     bitsByNumOn;
        std::vector<std::vector<size_t> >                   statesToBitsByNumOn;
        std::map< std::vector<size_t>, size_t >             bitsToStatesByNumOn;
        std::vector<std::string>                            stateDescriptionsByNumOn;
        std::vector<size_t>                                 rangeSizes;
        
        HomologousDiscreteCharacterData<StandardState>  data;
        std::string                                     format;
        size_t                                          num_taxa;
        size_t                                          numChars;
        size_t                                          numStates;
        size_t                                          numAllStates;
        
    };
}

#endif /* BitsetCharacterDataConverter_h */
