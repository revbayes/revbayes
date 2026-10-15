#ifndef TreeDistributionProperties_H
#define TreeDistributionProperties_H

namespace RevBayesCore {

    /**
     * @brief Properties of a tree distribution that a tree move can ask about.
     */
    class TreeDistributionProperties {

    public:
        virtual ~TreeDistributionProperties(void) {}

        virtual bool    allowsSampledAncestors(void) const { return false; }    //!< Checks if distribution is compatible with sampled ancestors
    };

}

#endif
