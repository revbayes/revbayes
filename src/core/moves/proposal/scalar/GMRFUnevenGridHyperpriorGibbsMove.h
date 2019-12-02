#ifndef GMRFUnevenGridHyperpriorGibbsMove_H
#define GMRFUnevenGridHyperpriorGibbsMove_H

#include "AbstractGibbsMove.h"
#include "AbstractMove.h"
#include "StochasticNode.h"

#include <set>
#include <vector>


namespace RevBayesCore {
    
    class GMRFUnevenGridHyperpriorGibbsMove : public AbstractGibbsMove {
        
    public:
        GMRFUnevenGridHyperpriorGibbsMove( StochasticNode<double> *g, std::vector< StochasticNode<double> *> n, TypedDagNode< RbVector<double> > *u, double z, double w);                                                         //!<  constructor
        virtual                                ~GMRFUnevenGridHyperpriorGibbsMove(void);                                                             //!< Destructor
        
        // Basic utility functions
        GMRFUnevenGridHyperpriorGibbsMove*               clone(void) const;                                                                  //!< Clone object
        const std::string&                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        
    protected:
        void                                    performGibbsMove(void);                                 //!< Perform move
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                             //!< Swap the pointers to the variable on which the move works on.
        
    private:
        
        // member variables
        StochasticNode<double>*                 global_scale;                                          //!< The half-Cauchy variable the Proposal is working on
        std::vector< StochasticNode<double> *>  normals;                                               //!< The normal random variables on which this is a hyperprior
        TypedDagNode< RbVector<double> > *      grid;                                                  //!< Local spacing of field
        double                                  zeta;                                                  //!< Controls total scale of field, x is halfCauchy(0,zeta)
    };
    
}


#endif
