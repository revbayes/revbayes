#ifndef Clade_H
#define Clade_H

#include <cstddef>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>
#include <boost/optional.hpp>

#include "RbBitSet.h"
#include "Taxon.h"

namespace RevBayesCore {
    
    /**
     * Object describing clades.
     *
     * A clade is simply a container of the taxon names.
     * Hence, this class just provides some convenience methods but could be considered as
     * a string-vector.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-03-10, version 1.0
     */
    class Clade  {
        
    public:
                                                    Clade(void);                                                            //! Default constructor: empty clade of age 0.0
                                                    Clade(const Taxon &t, const RbBitSet &b = RbBitSet() );                 //!< Default constructor with optional index
                                                    Clade(const std::set<Taxon> &n, const RbBitSet &b = RbBitSet() );       //!< Default constructor with optional index
                                                    Clade(const std::vector<Taxon> &n, const RbBitSet &b = RbBitSet() );    //!< Default constructor with optional index
                                                    Clade(const RbBitSet &b, const std::vector<Taxon> &n );    //!< Default constructor with optional index

        virtual                                    ~Clade(void) {}
        
        std::vector<Taxon>::const_iterator          begin(void) const;
        std::vector<Taxon>::iterator                begin(void);
        std::vector<Taxon>::const_iterator          end(void) const;
        std::vector<Taxon>::iterator                end(void);
        // overloaded operators
        bool                                        operator==(const Clade &t) const;
        bool                                        operator!=(const Clade &t) const;
        bool                                        operator<(const Clade &t) const;
        bool                                        operator<=(const Clade &t) const;
        bool                                        operator>(const Clade &t) const;
        bool                                        operator>=(const Clade &t) const;

        
        // Basic utility functions
        Clade*                                      clone(void) const;                                          //!< Clone object
        
        // public methods
        void                                        addTaxon(const Taxon &t);                                   //!< Add a taxon to our list.
        bool                                        conflicts(const Clade& c) const;                            //!< Is there a conflict between the clades, i.e. a partial overlap.
        double                                      getAge(void) const;                                         //!< Get the age of this clade.
        const RbBitSet&                             getBitRepresentation(void) const;                           //!< Get the clade as a bit representation.
        const std::string&                          getCladeName(void) const;                                   //!< Get the name of the clade.
        const std::set<Taxon>&                      getMrca(void) const;                                        //!< Get the mrca taxon.
        int                                         getNumberMissingTaxa(void) const;                           //!< Get the number of missing taxa.
        size_t                                      getNumberOfTaxa(void) const;                                //!< Get the number of taxa.
        const std::vector<Clade>&                   getOptionalConstraints(void) const;                         //!< Get optional clade constraints
        std::vector<Taxon>&                         getTaxa(void);                                              //!< Get the taxon names.
        const std::vector<Taxon>&                   getTaxa(void) const;                                        //!< Get the taxon names.
        const Taxon&                                getTaxon(size_t i) const;                                   //!< Get a single taxon name.
        const std::string&                          getTaxonName(size_t i) const;                               //!< Get a single taxon name.
        bool                                        hasOptionalConstraints(void) const;                         //!< Has optional clade constraints
        bool                                        isNegativeConstraint(void) const;                           //!< Get negative constraint flag.
        bool                                        isNestedWithin(const Clade& c) const;                       //!< Is the provided clade nested within me?
        bool                                        overlaps(const Clade& c) const;                             //!< Does the provided clade overlap with me?
        std::set<Taxon>                             intersection(const Clade&) const;                           //!< Get the taxa that both clades have in common.
        void                                        resetTaxonBitset(const std::map<std::string, size_t> map);
        void                                        setAge(double a);                                           //!< Set the age of the clade.
        void                                        setAges(const std::vector<Taxon>& taxa);                         //!< Set the aged of the taxa based on this set.
        void                                        setBitRepresentation(const RbBitSet &b);
        void                                        setCladeName(const std::string& n);                         //!< Set the name of the clade.
        void                                        setOptionalConstraints(const std::vector<Clade>& c);        //!< Set optional clade constraints.
        void                                        setMrca(const std::set<Taxon>&);                            //!< Set the mrca taxon, if applicable.
        void                                        setNumberMissingTaxa(int n);                                //!< Set the number of missing taxa in this clade.
        void                                        setTaxonAge(size_t i, double age);                          //!< Set a single taxon's age.
        void                                        setNegativeConstraint(bool);                                //!< Set clade to be a negative constraint
        size_t                                      size(void) const;                                           //!< Get the number of taxa.
        std::string                                 toString(void) const;                                       //!< Convert this value into a string.
        
        // public TopologyNode functions
        
    private: 
        
        // members
        double                                      age = 0.0;
        RbBitSet                                    bitset;
        std::string                                 clade_name;
        int                                         num_missing = 0;
        std::set<Taxon>                             mrca;
        std::vector<Taxon>                          taxa;
        bool                                        is_negative_constraint = false;
        boost::optional<std::vector<Clade>>         optional_constraints;
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Clade& x);                             //!< Overloaded output operator

    // Partial order
    bool cladeWithin(const Clade& c1, const Clade& c2);

    // Strict weak order
    bool cladeSmaller(const Clade& c1, const Clade& c2);
    bool cladeBefore(const Clade& c1, const Clade& c2);
}

#endif
