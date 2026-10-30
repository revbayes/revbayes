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
        TimeInterval(double min, double max);
        
        // Compare endpoints, ordered by maximum first and minimum second.
        bool                                operator==(const TimeInterval &t) const;
        std::partial_ordering               operator<=>(const TimeInterval &t) const;

        // Inclusive membership and point-interval detection.
        bool                                contains(double age) const;
        bool                                zeroWidth() const;

        // public functions
        double                              getMin(void) const;                                     //!< Get the lower endpoint.
        double                              getMax(void) const;                                     //!< Get the upper endpoint.
        
    private:
        // Explicit endpoints: neither is NaN. Infinity is allowed; min must not exceed max.
        double  min;
        double  max;
        
    };

    // Return the smallest interval containing both intervals.
    TimeInterval coveringInterval(const TimeInterval& a, const TimeInterval& b);

}

#endif
