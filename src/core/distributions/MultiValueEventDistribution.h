#ifndef MultiValueEventDistribution_H
#define MultiValueEventDistribution_H

#include <iosfwd>
#include <cstdint>
#include <vector>

#include "MultiValueEvent.h"
#include "TypedDistribution.h"
#include "MemberObject.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
    
    
    /**
     * This class implements a generic mixture distribution between several possible values.
     *
     * This mixture can be considered as a multinomial distribution. We specify a vector of probabilities
     * and a vector of values. Then, a value drawn from this distribution takes each value corresponding to
     * its probability.
     * The values are already of the correct mixture type. You may want to apply a mixture allocation move
     * to change between the current value. The values themselves change automatically when the input parameters change.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-11-18, version 1.0
     */
    class MultiValueEventDistribution : public TypedDistribution< MultiValueEvent >, public MemberObject< RbVector<double> >, public MemberObject< std::int64_t > {
        
    public:
        // constructor(s)
        MultiValueEventDistribution(TypedDistribution<std::int64_t> *ep, const std::vector< TypedDistribution<double> *> &vp, const std::vector< std::string > &n, const std::vector< std::int64_t > &min);
        MultiValueEventDistribution(const MultiValueEventDistribution &d);
        
        virtual                                            ~MultiValueEventDistribution();
        
        MultiValueEventDistribution&                        operator=(const MultiValueEventDistribution &d);

        // public member functions
        MultiValueEventDistribution*                        clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, std::int64_t &rv) const;     //!< Map the member methods to internal function calls
        const std::vector< std::int64_t >&                          getMinimumNumberOfEvents(void) const;
        const std::vector< TypedDistribution<double>* >&    getValuePriors(void) const;
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
        
    };
    
}



#endif

