#ifndef AutocorrelatedEventDistribution_H
#define AutocorrelatedEventDistribution_H

#include <iosfwd>
#include <vector>

#include "MultiValueEvent.h"
#include "TypedDistribution.h"
#include "MemberObject.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
template <class valueType> class RbVector;
    
    
    /**
     * This class implements a very generic approach to specify a process where the number of events is drawn from a base distribution.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    class AutocorrelatedEventDistribution : public TypedDistribution< MultiValueEvent >, public MemberObject< RbVector<double> >, public MemberObject< std::int64_t > {
        
    public:
        
        enum Autocorrelation { NONE, ACN, ACLN };
        
        // constructor(s)
        AutocorrelatedEventDistribution(TypedDistribution<std::int64_t> *ep, const std::vector< TypedDistribution<double> *>& vp, const std::vector< Autocorrelation >& ac, const std::vector< std::string >& ac_dep_var, const TypedDagNode< RbVector<double> >* ac_sd, const std::vector< std::string >& n, const std::vector< std::int64_t >& min, const std::string& sort_var);
        AutocorrelatedEventDistribution(const AutocorrelatedEventDistribution &d);
        
        virtual                                            ~AutocorrelatedEventDistribution();
        
        AutocorrelatedEventDistribution&                    operator=(const AutocorrelatedEventDistribution &d);

        // public member functions
        AutocorrelatedEventDistribution*                    clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;     //!< Map the member methods to internal function calls
        const std::vector< std::int64_t >&                          getMinimumNumberOfEvents(void) const;
        const std::vector< TypedDistribution<double>* >&    getValuePriors(void) const;
        bool                                                isAutocorrelated(size_t i) const;
        bool                                                isSorted(size_t i) const;
        void                                                redrawValue(void);
//        void                                                setNumberOfEvents(std::int64_t n);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        
        // helper methods
        void                                                simulate();
        
        // private members
        TypedDistribution<std::int64_t>*                            event_prior;
        std::vector< std::int64_t >                                 min_events;
        std::vector< std::string >                          names;
        std::vector< TypedDistribution<double>* >           value_priors;
        std::vector< Autocorrelation >                      autocorrelation_types;
        std::string                                         name_of_var_to_sort_by;
        int                                                 index_of_var_to_sort_by;
        std::vector< int >                                  autocorrelation_time_indeces;
        const TypedDagNode< RbVector<double> >*             autocorrelation_sigmas;

    };
    
}



#endif
