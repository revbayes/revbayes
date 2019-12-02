#ifndef DirichletDistribution_H
#define DirichletDistribution_H

#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Dirichlet distribution class.
     *
     * The Dirichlet distribution represents a family of distributions
     * on simplices (a vector that sums to one). The Dirichlet distribution has 1 parameters:
     *   alpha .. the relative counts in each category
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class DirichletDistribution : public TypedDistribution< Simplex > {
        
    public:
        DirichletDistribution(const TypedDagNode< RbVector<double> > *l);
        virtual                                            ~DirichletDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        DirichletDistribution*                              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode< RbVector<double> >*             alpha;
    };
    
}

#endif
