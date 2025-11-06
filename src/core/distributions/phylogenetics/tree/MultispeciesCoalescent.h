#ifndef MultispeciesCoalescent_H
#define MultispeciesCoalescent_H

#include "AbstractMultispeciesCoalescent.h"

namespace RevBayesCore {

    class Clade;

    class MultispeciesCoalescent : public AbstractMultispeciesCoalescent {

    public:
        MultispeciesCoalescent(const TypedDagNode<Tree> *st, const std::vector<Taxon> &t);
        virtual                                            ~MultispeciesCoalescent(void);                                                                       //!< Virtual destructor

        // public member functions
        MultispeciesCoalescent*                             clone(void) const;                                                                                  //!< Create an independent clone
        void                                                setNes(TypedDagNode<RbVector<double> >* inputNes);
        void                                                setNe(TypedDagNode<double>* inputNe);


    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        double                                              computeLnCoalescentProbability(size_t k, const std::vector<double> &t, double a, double b, size_t index, bool f);
        double                                              drawNe(size_t index);
        double                                              getNe(size_t index) const;


    private:

        // members
        const TypedDagNode<double >*                        Ne;
        const TypedDagNode<RbVector<double> >*              Nes;

    };

}

#endif
