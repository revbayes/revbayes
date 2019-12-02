#ifndef Split_H
#define Split_H

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Clade.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "Taxon.h"

namespace RevBayesCore {

    /**
     * Object describing subsplits.
     *
     * This is a bit-vector representation of a split in an unrooted tree.
     *
     * @copyright Copyright 2009-
     * @author Andrew Magee
     * @since 2018-11-23, version 1.0
     */
    class Split  {

    public:
                                                    Split(void);                                                            //! Default constructor: empty Split
                                                    Split(const std::vector<Taxon> &c, const std::vector<Taxon> &n );      //!< constructor: real Split
                                                    Split(RbBitSet &c);                               //!< Constructor from bitset representation of a clade

        virtual                                    ~Split(void) {}

        // overloaded operators
        Split&                                   operator=(const Split &t);
        bool                                        operator==(const Split &t) const;
        bool                                        operator!=(const Split &t) const;
        bool                                        operator<(const Split &t) const;
        bool                                        operator<=(const Split &t) const;
        bool                                        operator>(const Split &t) const;
        bool                                        operator>=(const Split &t) const;


        // Basic utility functions
        Split*                                      clone(void) const;                                          //!< Clone object

        // public methods
        RbBitSet                                    getBitset(void) const;                                      //!< Get the paired bitset representation of this subsplit as a clade
        bool                                        isCompatible(const Split &s) const;                      //!< Is argument subsplit compatible with this one?
        size_t                                      size(void) const;                                           //!< Get the number of taxa.
        std::string                                 toString(void) const;                                       //!< Convert this value into a string.

        // public TopologyNode functions

    private:

        // members
        RbBitSet                                    bitset;
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Split& x);                             //!< Overloaded output operator

}

#endif
