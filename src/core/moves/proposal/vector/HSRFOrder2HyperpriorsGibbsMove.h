#ifndef HSRFOrder2HyperpriorsGibbsMove_H
#define HSRFOrder2HyperpriorsGibbsMove_H

#include "AbstractGibbsMove.h"
#include "AbstractMove.h"
#include "StochasticNode.h"

#include <set>
#include <vector>

namespace RevBayesCore {

    class HSRFOrder2HyperpriorsGibbsMove : public AbstractGibbsMove {

    public:
        HSRFOrder2HyperpriorsGibbsMove( StochasticNode<double> *g, std::vector< StochasticNode<double> *> l, std::vector< StochasticNode<double> *> n, double z, double w);                                                         //!<  constructor
        virtual                                ~HSRFOrder2HyperpriorsGibbsMove(void);                                                             //!< Destructor

        // Basic utility functions
        HSRFOrder2HyperpriorsGibbsMove*            clone(void) const;                                                                  //!< Clone object
        const std::string&                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        bool                                     heatsAreAllowable(double prHeat, double lHeat, double pHeat); //<! A check we can call for whether it is okay to use this move with given prior/likelihood/posterior heats

    protected:
        void                                    performGibbsMove(void);                                 //!< Perform move
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                             //!< Swap the pointers to the variable on which the move works on.

    private:

        // member variables
        StochasticNode<double>*                 global_scale;                                          //!< The half-Cauchy variable the Proposal is working on
        std::vector< StochasticNode<double> *>  local_scales;                                          //!< The half-Cauchy variables the Proposal is working on
        std::vector< StochasticNode<double> *>  normals;                                               //!< The normal random variables on which this is a hyperprior
        double                                  zeta;                                                  //!< Controls total scale of field, x is halfCauchy(0,zeta)

    };

}

#endif
