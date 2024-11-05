#ifndef MultiValueEvent_h
#define MultiValueEvent_h

#include <cstdio>
#include <vector>
#include <iosfwd>

#include "RbVector.h"

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
    class MultiValueEvent {
        
    public:
        
        MultiValueEvent(void);                                              //!< Default constructor
        virtual                             ~MultiValueEvent() {}
        
        bool                                operator==(const MultiValueEvent &mve) const;
        bool                                operator!=(const MultiValueEvent &mve) const;
        bool                                operator<(const MultiValueEvent &mve) const;
        bool                                operator<=(const MultiValueEvent &mve) const { return operator<(mve) || operator==(mve); }

        // public methods
        void                                addValues(const RbVector<double> &v, const std::string &n);                               //!< Set the age.
        void                                clear(void);
        MultiValueEvent*                    clone(void) const;
        const std::string&                  getName(size_t i) const;
        std::int64_t                                getNumberOfEvents(void) const;
        size_t                              getNumberOfValues(void) const;
        RbVector<double>&                   getValues(size_t i);                                                //!< Get the values for this element.
        const RbVector<double>&             getValues(size_t i) const;                                          //!< Get the values for this element.
        RbVector<double>&                   getValues(const std::string& n);                              //!< Get the values for this element.
        const RbVector<double>&             getValues(const std::string& n) const;                              //!< Get the values for this element.
        void                                setNumberOfEvents(std::int64_t n);
        void                                setValues(const RbVector<double>& v, const std::string& n);                               //!< Set the age.


    private:
        
        // private members
        std::vector<std::string>            names;
        std::int64_t                                num_events;
        std::vector<RbVector<double> >      values;

    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MultiValueEvent& x);                                         //!< Overloaded output operator
}


#endif /* MultiValueEvent_hpp */
