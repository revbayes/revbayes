#ifndef MultispeciesCoalescentInverseGammaPrior_H
#define MultispeciesCoalescentInverseGammaPrior_H

#include "AbstractMultispeciesCoalescentGenewise.h"

namespace RevBayesCore {

    class Clade;

    class MultispeciesCoalescentInverseGammaPrior : public AbstractMultispeciesCoalescentGenewise {

    public:
        MultispeciesCoalescentInverseGammaPrior(const TypedDagNode<Tree> *st, RbVector< RbVector<Taxon> > t, size_t ngt);
        virtual                                            ~MultispeciesCoalescentInverseGammaPrior(void);                                                                       //!< Virtual destructor

        // public member functions
        MultispeciesCoalescentInverseGammaPrior*            clone(void) const;                                                                                  //!< Create an independent clone
        void                                                setShape(TypedDagNode<double>* s);
        void                                                setRate(TypedDagNode<double>* r);


    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        double                                              computeLnCoalescentProbability(std::vector<size_t> k, const std::vector< std::vector<double> > &t, double a, double b, size_t index, bool f);
        double                                              drawNe(size_t index);


    private:

        // members
        const TypedDagNode<double>*                          shape;
        const TypedDagNode<double>*                          rate;

    };

}

#endif
