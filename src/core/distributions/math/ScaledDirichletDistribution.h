#ifndef ScaledDirichletDistribution_H
#define ScaledDirichletDistribution_H

#include "Simplex.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief ScaledDirichlet distribution class.
     *
     * The ScaledDirichlet distribution represents a family of distributions
     * on simplices (a vector that sums to one). The ScaledDirichlet distribution has 2 parameters:
     *   alpha .. the relative counts in each category
     * Instances of this class can be associated to stochastic variables.
     *
     *
     */
    class ScaledDirichletDistribution : public TypedDistribution<Simplex> {
        
    public:
        ScaledDirichletDistribution(const TypedDagNode< RbVector<double> > *a, const TypedDagNode< RbVector<double> > *b );
        virtual                                            ~ScaledDirichletDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        ScaledDirichletDistribution*                        clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             alpha;
        const TypedDagNode< RbVector<double> >*             beta;
    };
    
}

#endif
