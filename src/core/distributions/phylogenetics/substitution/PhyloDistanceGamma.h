#ifndef PhyloDistanceGamma_H
#define PhyloDistanceGamma_H

#include "DistanceMatrix.h"
#include "MatrixReal.h"
#include "TopologyNode.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class PhyloDistanceGamma : public TypedDistribution< DistanceMatrix >/*, public TreeChangeEventListener*/ {
        
    public:
		PhyloDistanceGamma( const TypedDagNode< Tree > *t );//, std::vector<std::string> names ); //, MatrixReal distMatrix, MatrixReal varMatrix );
		
        virtual                                            ~PhyloDistanceGamma(void);                                                                   //!< Virtual destructor
        
        // public member functions
                PhyloDistanceGamma*                                                 clone(void) const;                                                                          //!< Create an independent clone
		double                                                              computeLnProbability(void);
		//void                                                                fireTreeChangeEvent(const TopologyNode &n);                                                 //!< The tree has changed and we want to know which part.
		void                                                                redrawValue(void);
		void                                                                setDistanceMatrix(const TypedDagNode< DistanceMatrix > *dm);
		void                                                                setVarianceMatrix(const TypedDagNode< DistanceMatrix > *dm);
		void                                                                setNames(const std::vector< std::string >& n);
		//void 																setValue(const MatrixReal *dm, const MatrixReal *vm); //!< Set the current value, e.g. attach an observation (clamp)
		void 																simulate( );


		
    protected:
		
		// Parameter management functions.
		virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                                    //!< Swap a parameter
		
		void 																updateAlphaAndBetaMatrices();
		
		// members
		const TypedDagNode<Tree>*                                           tau;
		size_t 																numTips;
		const TypedDagNode<DistanceMatrix>*                                 distanceMatrix;
		const TypedDagNode<DistanceMatrix>*                                 varianceMatrix;
		DistanceMatrix                                                      alphaMatrix;
		DistanceMatrix                                                      betaMatrix;
		std::vector<std::string> 											matrixNames;
		double 																lnProb;

		// the likelihoods
		std::vector<double>                                                 partialLikelihoods;
		

    private:
		double 																computeLogLikelihood( void );

		
    
    };
    
}


#endif
