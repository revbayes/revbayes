#ifndef DPPGibbsConcentrationMove_H
#define DPPGibbsConcentrationMove_H

#include <iosfwd>

#include "AbstractGibbsMove.h"


namespace RevBayesCore {
class DagNode;
template <class valueType> class DeterministicNode;
template <class valueType> class TypedDagNode;
template <class variableType> class StochasticNode;

    /**
     * @brief The Move for updating the concentration parameter of the DPP
     *
     * This is a Gibbs move on the concentration parameter (alpha) where the prior on
     * alpha must be a Gamma distribution with fixed hyperparameters.
     * This move is described by Escobar & West (Bayesian density estimation and inference using mixtures, 1995),
     * (also see Dozario, On selecting a prior for the precision parameter of Dirichlet process mixture models, 2009).
     * This move is conditional on the # of categories in the current state of the DPP
     *
     * The move has the following member variables:
     *  @param variable the concentration parameter (double)
     *  @param numCats the number of mixture categories in the current state of the DPP (std::int64_t)
     *  @param gammaShape the shape parameter of the gamma prior (double)
     *  @param gammaRate the rate parameter of the gamma prior (double)
     *  @param numElem the number of elements in the DPP (int)
     *
     */

    class DPPGibbsConcentrationMove : public AbstractGibbsMove {
		
    public:
        DPPGibbsConcentrationMove(StochasticNode<double> *n, DeterministicNode<std::int64_t>* v, TypedDagNode< double >* gS, TypedDagNode< double >* gR, int ne, double w);                                                                      //!< Internal constructor
		
        // Basic utility functions
        DPPGibbsConcentrationMove*                              clone(void) const;                                                                  //!< Clone object
        const std::string&                                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
		
    protected:
        void													performGibbsMove(void);                                                            //!< Perform move
        void                                                    swapNodeInternal(DagNode *oldN, DagNode *newN);
        
    private:
		
        StochasticNode<double>*									variable;
		DeterministicNode<std::int64_t>*								numCats;
		TypedDagNode< double >*									gammaShape;
		TypedDagNode< double >*									gammaRate;
		int														numElem;
    };
}



#endif

