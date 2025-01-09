#ifndef Taxon_H
#define Taxon_H

#include <ostream>
#include <map>

#include "TimeInterval.h"

namespace RevBayesCore {
    
    
    /**
     * @brief Taxon object representing an individual.
     *
     * This class represents a taxon.
     * A taxon simply consists of a taxon name and a list of occurrences when this taxon has been observed.
     * Occurrences are stored in a map associating an age range with an occurrence count.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     */
    class Taxon {
        
    public:
        
                                              Taxon(void);                                //!< Default constructor
                                              Taxon(const std::string &n);                //!< Regular constructor
        virtual                               ~Taxon(void) {}
        
        bool                                  operator==(const Taxon &t) const;           //!< Equals operators
        bool                                  operator!=(const Taxon &t) const;           //!< Not-quals operators
        bool                                  operator<(const Taxon &t) const;            //!< Less-than operators
        bool                                  operator<=(const Taxon &t) const;           //!< Less-than operators
        bool                                  operator>(const Taxon &t) const;            //!< Less-than operators
        bool                                  operator>=(const Taxon &t) const;           //!< Less-than operators
        
        // public methods
        void                                  addOccurrence(const TimeInterval &d);       //!< Add an occurrence.
        bool                                  isExtinct(void) const;                      //!< Get the extinct status flag.
        double                                getAge(void) const;                         //!< Get the age for this taxon.
        const TimeInterval&                   getAgeRange(void) const;                    //!< Get the date info for this taxon.
        const std::string                     getJsonRespresentation(void) const;         //!< Get JSON-formatted string.
        double                                getMaxAge(void) const;                      //!< Get the max age for this taxon.
        double                                getMinAge(void) const;                      //!< Get the min age for this taxon.
        const std::string&                    getName(void) const;                        //!< Get the name for this taxon.
        const std::map<TimeInterval, size_t>& getOccurrences(void) const;                 //!< Get the occurrence ages and counts.
        const std::string&                    getSpeciesName(void) const;                 //!< Get the name of the species.

        void                                  setAge(double a);                           //!< Set the age.
        void                                  setMinAge(double a);                        //!< Set the min age.
        void                                  setMaxAge(double a);                        //!< Set the max age.
        void                                  setAgeRange(const TimeInterval &d);         //!< Set the date info.
        void                                  setExtinct(bool extinct);                   //!< Set the extinct status flag.
        void                                  setName(const std::string &n);              //!< Set the name.
        void                                  setSpeciesName(const std::string &n);       //!< Set the name of the species.
        
    private:
        
        // private members
        TimeInterval                          age_range;
        std::map<TimeInterval, size_t>        occurrences;
        std::string                           name;
        std::string                           species_name;
        bool                                  extinct;
    
    };

    // Global functions using the class
    std::ostream&                             operator<<(std::ostream& o, const Taxon& x);//!< Overloaded output operator
}

#endif
