#ifndef PhyloCTMCSiteHomogeneousDolloBinary_H
#define PhyloCTMCSiteHomogeneousDolloBinary_H

#include <cstddef>
#include <map>
#include <vector>

#include "PhyloCTMCSiteHomogeneousDollo.h"

namespace RevBayesCore {
class Tree;
template <class valueType> class TypedDagNode;

    struct DolloBinaryAscertainmentBias {

      enum Coding { NOABSENCESITES      = 0x01,
                    NOPRESENCESITES     = 0x02,
                    VARIABLE            = 0x03,
                    NOSINGLETONPRESENCE = 0x04,
                    NOSINGLETONABSENCE  = 0x08,
                    NOSINGLETONS        = 0x0C,
                    INFORMATIVE         = 0x0F
                  };
    };

    class PhyloCTMCSiteHomogeneousDolloBinary : public PhyloCTMCSiteHomogeneousDollo {

        public:
        PhyloCTMCSiteHomogeneousDolloBinary(const TypedDagNode< Tree > *t, bool c, size_t nSites, bool amb, DolloAscertainmentBias::Coding cod = DolloAscertainmentBias::NOABSENCESITES, bool normalize = true);

        protected:
            bool                                        isSitePatternCompatible( std::map<size_t, size_t> );
            double                                      computeIntegratedNodeCorrection(const std::vector<std::vector<std::vector<double> > >& partials, size_t nodeIndex, size_t mask, size_t mixture, const std::vector<double> &f);
        };

}

#endif
