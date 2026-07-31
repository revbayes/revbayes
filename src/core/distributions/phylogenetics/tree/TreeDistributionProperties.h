#ifndef TreeDistributionProperties_H
#define TreeDistributionProperties_H

#include "Taxon.h"

namespace RevBayesCore {

    /**
     * @brief What a tree move or a tree distribution's own machinery may assume about the process.
     *
     * A move holds only a TypedDistribution<Tree>, and the two tree distributions that answer these
     * are siblings rather than one deriving from the other: AbstractRootedTreeDistribution, and
     * TopologyConstrainedTreeDistribution, which wraps one. This is the interface they share, so a
     * question can be asked without it living on the distribution base for every value type, and
     * without a wrapped distribution going unheard.
     *
     * isExtended is the property; the other two are policies that default to it and exist so a
     * process with a finer answer can give one. FossilizedBirthDeathSpeciationProcess gives all
     * three, since which of its tips are pinned to a fossil range is a per-taxon question.
     */
    class TreeDistributionProperties {

    public:
        virtual ~TreeDistributionProperties(void) {}

        virtual bool    isExtended(void) const { return false; }                                     //!< Extended trees have tips at extinctions, which may sit below a fossil's age range
        virtual bool    allowsSampledAncestors(void) const { return false; }                         //!< Whether a sampled ancestor is a configuration this distribution can score
        virtual bool    tipAgeConstrainedToRange(const Taxon &t) const { return isExtended() == false; }
        virtual bool    validatesTipAgesOnSet(void) const { return isExtended() == false; }
    };

}

#endif
