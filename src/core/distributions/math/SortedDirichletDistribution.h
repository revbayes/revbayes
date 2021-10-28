#ifndef SortedDirichletDistribution_H
#define SortedDirichletDistribution_H

#include "Simplex.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;

/**
* @brief Sorted Dirichlet distribution class.
*
* This class is an adaptation to the Dirichlet distribution that has the additional constraint that all values
* in the distribution are sorted in descending order.
*
*/

    class SortedDirichletDistribution : public TypedDistribution< Simplex > {

        public:
        SortedDirichletDistribution(const TypedDagNode< RbVector<double> > *l);  //!< Constructor
        virtual ~SortedDirichletDistribution(void);                              //!< Virtual destructor

        SortedDirichletDistribution* clone(void) const;                          //!< Create an independent clone
        double                       computeLnProbability(void);                 //!< Natural log of the probability density
        void                         redrawValue(void);
        
    protected:
        void                         swapParameterInternal(const DagNode *oldP, const DagNode *newP); //!< Swap a parameter
        
    private:
        const TypedDagNode< RbVector<double> >*  alpha;
    };
}

#endif
