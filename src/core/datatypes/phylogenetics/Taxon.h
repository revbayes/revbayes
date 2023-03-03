#ifndef Taxon_H
#define Taxon_H

#include <ostream>
#include <boost/optional.hpp>

#include "TimeInterval.h"

namespace RevBayesCore {
    
    
    /**
     * @brief Taxon object representing an individual.
     *
     * This class represents a taxon.
     * A taxon simply consists of a taxon name and a date when this taxon has been observed.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     */
    class Taxon {
        
    public:
        
                                            Taxon(void);                                //!< Default constructor
                                            Taxon(const std::string &n);                //!< Regular constructor
        virtual                             ~Taxon(void) {}
        
        bool                                operator==(const Taxon &t) const;           //!< Equals operators
        bool                                operator!=(const Taxon &t) const;           //!< Not-quals operators
        bool                                operator<(const Taxon &t) const;            //!< Less-than operators
        bool                                operator<=(const Taxon &t) const;           //!< Less-than operators
        bool                                operator>(const Taxon &t) const;            //!< Less-than operators
        bool                                operator>=(const Taxon &t) const;           //!< Less-than operators
        
        // public methods
        double                              getAge(void) const;                         //!< Get the age for this taxon.
        const TimeInterval&                 getAgeRange(void) const;                    //!< Get the date info for this taxon.
        const std::string                   getJsonRespresentation(void) const;         //!< Get JSON-formatted string
        const std::string&                  getName(void) const;                        //!< Get the name for this taxon.
        const std::string&                  getPloidy(void) const;                      //!< Get the ploidy level for this taxon.
        const std::string&                  getSpeciesName(void) const;                 //!< Get the name of the species.
        void                                setAge(double a);                           //!< Set the age.
        void                                setAgeRange(const TimeInterval &d);         //!< Set the date info.
        void                                setName(const std::string &n);              //!< Set the name.
        void                                setPloidy(const std::string &n);            //!< Set the ploidy information.
        void                                setSpeciesName(const std::string &n);       //!< Set the name of the species.
        
    private:
        
        // private members
        TimeInterval                        age_range;
        std::string                         name;
        std::string                         species_name;
        boost::optional<std::string>        ploidy;
    
    };
	
    // Global functions using the class
    std::ostream&                           operator<<(std::ostream& o, const Taxon& x);//!< Overloaded output operator
}

#endif
