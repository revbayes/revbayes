#ifndef TimeInterval_H
#define TimeInterval_H

#include <compare>

namespace RevBayesCore {
    
    /**
     * @brief Time interval object for dates of taxon samples.
     *
     * This class provides the time interval information for a taxon.
     * Some convenience methods are provided.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     */
    class TimeInterval {
        
    public:
        TimeInterval(); //!< Default constructor
        TimeInterval(double min, double max);
        
        // Compare effective endpoints, ordered by maximum first and minimum second.
        bool                                operator==(const TimeInterval &t) const;
        std::partial_ordering               operator<=>(const TimeInterval &t) const;

        // Membership and zero width use effective endpoints, including the current [0, 0] default.
        bool                                contains(double age) const;
        bool                                zeroWidth() const;

        // public functions
        double                              getMin(void) const;                                     //!< Get the number of days since counting.
        double                              getMax(void) const;                                     //!< Get the number of million years ago.
        void                                setMin(double);                                         //!< Get the number of days since counting.
        void                                setMax(double);                                         //!< Get the number of million years ago.
        
    private:
        double  min;
        double  max;
        
    };

    // Return the smallest interval containing both effective intervals, with explicit bounds.
    TimeInterval coveringInterval(const TimeInterval& a, const TimeInterval& b);

}

#endif
