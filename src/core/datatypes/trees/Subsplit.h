#ifndef Subsplit_H
#define Subsplit_H

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
     * A subsplit is simply a pair of clades.
     *
     * @copyright Copyright 2009-
     * @author Andrew Magee
     * @since 2018-11-23, version 1.0
     */
    class Subsplit  {

    public:
                                                    Subsplit(void);                                                            //! Default constructor: empty Subsplit
                                                    Subsplit(const std::vector<Taxon> &c1, const std::vector<Taxon> &n );                                                //!< Default constructor: fake Subsplit
                                                    Subsplit(const std::vector<Taxon> &c1, const std::vector<Taxon> &c2, const std::vector<Taxon> &n );                               //!< Default constructor: real Subsplit
                                                    Subsplit(const RbBitSet &c1, const RbBitSet &c2 );                               //!< Default constructor: real Subsplit from bitsets
                                                    Subsplit(const Subsplit &s1, const Subsplit &s2 );                               //!< Default constructor: Subsplit from its child subsplits

        virtual                                    ~Subsplit(void) {}

        // overloaded operators
        Subsplit&                                   operator=(const Subsplit &t);
        bool                                        operator==(const Subsplit &t) const;
        bool                                        operator!=(const Subsplit &t) const;
        bool                                        operator<(const Subsplit &t) const;
        bool                                        operator<=(const Subsplit &t) const;
        bool                                        operator>(const Subsplit &t) const;
        bool                                        operator>=(const Subsplit &t) const;


        // Basic utility functions
        Subsplit*                                      clone(void) const;                                          //!< Clone object

        // public methods
        Clade                                       asClade(const std::vector<Taxon> &taxa) const;                                        //!< Get clades Y and Z as a single clade
        RbBitSet                                    asCladeBitset(void) const;                                  //!< Convert this value into a bitset representation for a clade.
        RbBitSet                                    asSplitBitset(void) const;                                        //!< Convert this value into a split (as bitset representation).
        std::vector<std::pair<Subsplit,Subsplit> >  doVirtualRootingNonRootParent(const Subsplit &parent, const Subsplit &child) const; //!< Gives us all parent-child subsplits that will replace this one on the equivalent edge in a tree rerooted to a specific edge
        void                                        doVirtualRootingNonRootParent(const Subsplit &parent, const Subsplit &child, std::vector<std::pair<Subsplit,Subsplit> > &cases) const; //!< Gives us all parent-child subsplits that will replace this one on the equivalent edge in a tree rerooted to a specific edge
        std::vector<std::pair<Subsplit,Subsplit> >  doVirtualRootingRootParent(const Subsplit &sister1, const Subsplit &sister2, const Subsplit &child) const; //!< Gives us all parent-child subsplits that will replace this one on the equivalent edge in a tree rerooted to a specific edge
        void                                        doVirtualRootingRootParent(const Subsplit &sister1, const Subsplit &sister2, const Subsplit &child, std::vector<std::pair<Subsplit,Subsplit> > &cases) const; //!< Gives us all parent-child subsplits that will replace this one on the equivalent edge in a tree rerooted to a specific edge
        std::pair<RbBitSet,RbBitSet>                getBitset(void) const;                                      //!< Get the paired bitset representation of this subsplit as a bitset
        const std::pair<RbBitSet,RbBitSet>&         getBitsetRef(void) const;                                      //!< Get the paired bitset representation of this subsplit as a bitset referemce
        const std::pair<RbBitSet,RbBitSet>*         getBitsetPtr(void) const;                                      //!< Get the paired bitset representation of this subsplit as a bitset referemce
        size_t                                      getFsbY(void) const;                                        //!< Get first set bit in Y
        size_t                                      getFsbZ(void) const;                                        //!< Get first set bit in Z
        Clade                                       getY(const std::vector<Taxon>& taxa) const;                                           //!< Get clade Y
        RbBitSet                                    getYBitset(void) const;                                     //!< Get clade Y as bitset
        RbBitSet&                                   getYBitsetRef(void) const;                                     //!< Get clade Y as bitset  reference
        Clade                                       getZ(const std::vector<Taxon>& taxa) const;                                           //!< Get clade Z
        RbBitSet                                    getZBitset(void) const;                                     //!< Get clade Z as bitset
        RbBitSet&                                   getZBitsetRef(void) const;                                     //!< Get clade Z as bitset reference
        size_t                                      getHash(void) const;                                         //!< A hash function (mainly for unordered_maps)
        bool                                        isChildOfY(const Subsplit &s) const;                      //!< Is argument subsplit compatible with this one?
        bool                                        isChildOfZ(const Subsplit &s) const;                      //!< Is argument subsplit compatible with this one?
        bool                                        isCompatible(const Subsplit &s) const;                      //!< Is argument subsplit compatible with this one?
        bool                                        isFake() const;                                             //!< Is this a fake subsplit (aka a single tip)?
        Subsplit                                    rootSplitFromClade() const;                                 //!< If this clade (that this subsplit represents) is one side of a root split, make the full root subsplit.
        size_t                                      size(void) const;                                           //!< Get the number of taxa.
        bool                                        splitsAreDisjoint(void) const;
        std::string                                 toString(void) const;                                       //!< Convert this value into a string.
        std::pair<Subsplit,Subsplit>                reverseParentChildNonRootParent(const Subsplit &t, const Subsplit &s) const;       //!< Gives us the parent-child subsplit that will replace this one on the equivalent edge in a tree rerooted in the child
        std::pair<Subsplit,Subsplit>                reverseParentChildRootParent(const Subsplit &sister1, const Subsplit &sister2, const Subsplit &s) const;       //!< Gives us the parent-child subsplit that will replace this one on the equivalent edge in a tree rerooted in the child

        // public TopologyNode functions

    private:

        // members
        std::pair<RbBitSet,RbBitSet>                bitset;
        bool                                        is_fake;
        size_t                                      fsb_y;
        size_t                                      fsb_z;
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Subsplit& x);                             //!< Overloaded output operator

}

#endif
