#ifndef MultispeciesCoalescentUniformPrior_H
#define MultispeciesCoalescentUniformPrior_H

#include "AbstractMultispeciesCoalescent.h"

namespace RevBayesCore {
    
    class Clade;
    
    class MultispeciesCoalescentUniformPrior : public AbstractMultispeciesCoalescent {
        
    public:
        MultispeciesCoalescentUniformPrior(const TypedDagNode<Tree> *st, const std::vector<Taxon> &t);
        virtual                                            ~MultispeciesCoalescentUniformPrior(void);                                                                       //!< Virtual destructor
        
        // public member functions
        MultispeciesCoalescentUniformPrior*                 clone(void) const;                                                                                  //!< Create an independent clone
        void                                                setMaxTheta(TypedDagNode<double>* m);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        double                                              computeLnCoalescentProbability(size_t k, const std::vector<double> &t, double a, double b, size_t index, bool f);
        double                                              drawNe(size_t index);

        
        
    private:
        
        // members
        const TypedDagNode<double >*                        max_theta;
        
    };
    
}

#endif
