#ifndef PhyloCTMCSiteHomogeneousBinary_H
#define PhyloCTMCSiteHomogeneousBinary_H

#include <cstddef>
#include <map>
#include <vector>

#include "BinaryState.h"
#include "PhyloCTMCSiteHomogeneousConditional.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DiscreteTaxonData.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "RbException.h"
#include "TypedDistribution.h"
#include "RbBitSet.h"

namespace RevBayesCore {
class Tree;
template <class valueType> class TypedDagNode;

    struct BinaryAscertainmentBias {
        
      enum Coding { ALL                 = 0x00,
                    NOABSENCESITES      = 0x01,
                    NOPRESENCESITES     = 0x02,
                    VARIABLE            = 0x03,
                    NOSINGLETONPRESENCE = 0x04,
                    NOSINGLETONABSENCE  = 0x08,
                    NOSINGLETONS        = 0x0C,
                    INFORMATIVE         = 0x0F
                  };
    };

    class PhyloCTMCSiteHomogeneousBinary : public PhyloCTMCSiteHomogeneousConditional<BinaryState> {

        public:
            PhyloCTMCSiteHomogeneousBinary(const TypedDagNode< Tree > *t, bool c, size_t nSites, bool amb, BinaryAscertainmentBias::Coding cod = BinaryAscertainmentBias::ALL, bool internal = false, bool gapmatch = true, bool invert = false);

            // public member functions
            PhyloCTMCSiteHomogeneousBinary*                     clone(void) const;

        protected:

            double                                              sumRootLikelihood( void );
            bool                                                isSitePatternCompatible(std::map<size_t, size_t>);
            bool                                                isSitePatternCompatible(std::map<RbBitSet, size_t>);
        };

}

#endif
