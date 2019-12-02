/**
 * @file
 * This file contains the declaration of DPPGibbsConcentrationMove, which performs a fully conditional Gibbs move on the concentration parameter of the DPP
 *   this move is conditional on the # of categories in the current state of the DPP
 *
 * @brief Declaration of DPPGibbsConcentrationMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-11 14:54:35 +0200 (Fri, 11 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-20, version 1.0
 *
 * $Id: DPPGibbsConcentrationMove.h $
 */

#ifndef DPPGibbsConcentrationMove_H
#define DPPGibbsConcentrationMove_H

#include <ostream>

#include "AbstractGibbsMove.h"
#include "DirichletProcessPriorDistribution.h"
#include "StochasticNode.h"
#include "DeterministicNode.h"
#include "TypedDagNode.h"
#include "ConstantNode.h"


namespace RevBayesCore {
	
    class DPPGibbsConcentrationMove : public AbstractGibbsMove {
		
    public:
        DPPGibbsConcentrationMove(StochasticNode<double> *n, DeterministicNode<long>* v, TypedDagNode< double >* gS, TypedDagNode< double >* gR, int ne, double w);                                                                      //!< Internal constructor
		
        // Basic utility functions
        DPPGibbsConcentrationMove*                              clone(void) const;                                                                  //!< Clone object
        const std::string&                                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
		
    protected:
        void													performGibbsMove(void);                                                            //!< Perform move
        void                                                    swapNodeInternal(DagNode *oldN, DagNode *newN);
        
    private:
		
        StochasticNode<double>*									variable;
		DeterministicNode<long>*									numCats;
		TypedDagNode< double >*									gammaShape;
		TypedDagNode< double >*									gammaRate;
		int														numElem;
    };
}



#endif

