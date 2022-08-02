#ifndef TimeInterval_H
#define TimeInterval_H

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
        
        bool                                operator==(const TimeInterval &t) const;                //!< Equals operator
        bool                                operator!=(const TimeInterval &t) const;                //!< Not-equals operator
        bool                                operator<(const TimeInterval &t) const;                 //!< Less than operator

        // public functions
        double                              getMin(void) const;                                     //!< Get the number of days since counting.
        double                              getMax(void) const;                                     //!< Get the number of million years ago.
        void                                setMin(double);                                         //!< Get the number of days since counting.
        void                                setMax(double);                                         //!< Get the number of million years ago.
        
    private:
        double  min;
        double  max;
        
    };

}

#endif
