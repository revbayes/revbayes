#ifndef MultiValueEventDistribution_H
#define MultiValueEventDistribution_H

#include "MultiValueEvent.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    
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
    class MultiValueEventDistribution : public TypedDistribution< MultiValueEvent >, public MemberObject< RbVector<double> >, public MemberObject< long > {
        
    public:
        // constructor(s)
        MultiValueEventDistribution(TypedDistribution<long> *ep, const std::vector< TypedDistribution<double> *> &vp, const std::vector< std::string > &n, const std::vector< long > &min);
        MultiValueEventDistribution(const MultiValueEventDistribution &d);
        
        virtual                                            ~MultiValueEventDistribution();
        
        MultiValueEventDistribution&                        operator=(const MultiValueEventDistribution &d);

        // public member functions
        MultiValueEventDistribution*                        clone(void) const;                                                                                  //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, long &rv) const;     //!< Map the member methods to internal function calls
        const std::vector< long >&                          getMinimumNumberOfEvents(void) const;
        const std::vector< TypedDistribution<double>* >&    getValuePriors(void) const;
        void                                                redrawValue(void);
//        void                                                setNumberOfEvents(long n);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Swap a parameter
        
        
    private:
        
        // helper methods
        void                                                simulate();
        
        // private members
        TypedDistribution<long>*                            event_prior;
        std::vector< long >                                 min_events;
        std::vector< std::string >                          names;
        std::vector< TypedDistribution<double>* >           value_priors;
        
    };
    
}



#endif

