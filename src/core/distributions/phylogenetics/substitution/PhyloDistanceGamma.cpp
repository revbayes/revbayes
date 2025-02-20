#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "PhyloDistanceGamma.h"
#include "DistributionGamma.h"
#include "StochasticNode.h"
#include "DistanceMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


PhyloDistanceGamma::PhyloDistanceGamma( const TypedDagNode< Tree > *t ) : TypedDistribution<DistanceMatrix>( new DistanceMatrix( t->getValue().getNumberOfTips()) ),
    tau( t ),
    numTips(t->getValue().getNumberOfTips() ),
    distanceMatrix( new ConstantNode<DistanceMatrix>("distanceMatrix", new DistanceMatrix(numTips) ) ),
    varianceMatrix( new ConstantNode<DistanceMatrix>("varianceMatrix", new DistanceMatrix(numTips) ) ),
    alphaMatrix( numTips ),
    betaMatrix( numTips )
/*, distanceMatrix(distMatrix), varianceMatrix(varMatrix), matrixNames (names)*/
//    distanceMatrix(distMatrix), varianceMatrix(varMatrix), matrixNames (names)*/
{
    
    addParameter( tau );
    addParameter( distanceMatrix );
    addParameter( varianceMatrix );
    
    lnProb = 0.0 ;
    
    //fill alpha and beta matrices
    updateAlphaAndBetaMatrices( );
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloDistanceGamma::~PhyloDistanceGamma( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
  /*  if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }*/
    
}



PhyloDistanceGamma* PhyloDistanceGamma::clone( void ) const
{
    
    return new PhyloDistanceGamma( *this );
}


double PhyloDistanceGamma::computeLogLikelihood( void )
{
    
    //First, compute the distance matrix from the current tree
    //DistanceMatrix* mat  = RevBayesCore::TreeUtilities::getDistanceMatrix<treeType> ( tau->getValue() );
	
	std::cout << "computeLogLikelihood: " << (*(this->value) )[3][5] <<std::endl;
	delete this->value;
	this->value  = RevBayesCore::TreeUtilities::getDistanceMatrix ( tau->getValue() );
	
	std::cout << "computeLogLikelihood 2: " << (*(this->value) )[3][5] <<std::endl;

	
    // Now we need to know the order in which the distances have been put in the matrix.
    // We know they are output in the same order the tipnames are given.
    std::cout << "computeLogLikelihood "<<std::endl;
    std::vector< std::string > namesFromTree = tau->getValue().getTipNames();
    
    std::cout << "computeLogLikelihood 2"<<std::endl;
    //We build a map linking names to their id in the tree-based matrix.
    std::map<std::string, size_t > namesToId;
    for (size_t i = 0; i < namesFromTree.size() ; ++i) {
        namesToId [ namesFromTree[i] ] = i;
    }
    std::cout << "computeLogLikelihood 3"<<std::endl;
    
    //Second, for each pairwise distance, compute its logprobability according to a Gamma distribution with parameters alpha and beta.
    double alpha, beta;
    double logL = 0.0;
    std::string nameI;
    std::string nameJ;
    for (size_t i = 0; i < numTips - 1; ++i) {
        nameI =matrixNames[i];
        for (size_t j = i+1; j < numTips; ++j) {
            nameJ =matrixNames[j];
            alpha = alphaMatrix[i][j];
            beta = betaMatrix[i][j];
            logL += RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] );
			/*
			//various checks:
			 std::cout <<"alpha: "<< alpha << " beta: " << beta << "; (*(this->value) )[namesToId[nameI]][namesToId[nameJ]]: " << (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] << " ; lnPdf: " << RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] )<< std::endl;
			std::vector<double > meanV ;
			double mean = 0;
			double variance = 0;
			for (size_t n = 0; n<1000; n++) {
				meanV.push_back( RbStatistics::Gamma::rv( alpha, beta, *GLOBAL_RNG ) );
				mean += meanV[n];
				variance += meanV[n] * meanV[n];
			}
			mean = mean /1000;
			variance = variance/1000 - mean * mean;
			std::cout <<"alpha: "<< alpha << " beta: " << beta << "; mean: " << mean << "; variance: "<< variance <<std::endl;
			std::cout << "Attempts: " <<std::endl;
			std::cout << "-0.05: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] - 0.05 )<<std::endl;
			std::cout << "-0.02: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] - 0.02 )<<std::endl;
			std::cout << "-0.01: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] - 0.01 )<<std::endl;
			std::cout << "MAX ?: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]]  )<<std::endl;
			std::cout << "+0.01: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] + 0.01 )<<std::endl;
			std::cout << "+0.02: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] + 0.02 )<<std::endl;
			std::cout << "+0.05: "<< RbStatistics::Gamma::lnPdf( alpha, beta, (*(this->value) )[namesToId[nameI]][namesToId[nameJ]] + 0.05 )<<std::endl;
			*/
        }
    }
    //delete mat;
	std::cout << "computeLogLikelihood 4"<<std::endl;

    return logL;
}


double PhyloDistanceGamma::computeLnProbability( void )
{
    
    // sum the partials up
    this->lnProb = computeLogLikelihood();
    
    return this->lnProb;
    
}


/*
void PhyloDistanceGamma::fireTreeChangeEvent( const TopologyNode &n )
{
	
    // call a recursive flagging of all node above (closer to the root) and including this node
    // For the moment we do nothing...
}
*/


void PhyloDistanceGamma::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;
    
    // create a new MatrixReal object
    this->value = distanceMatrix->getValue().clone();
    //TreePairwiseDistances ( tau );
    
    
}

/*
 template<class Tree>
 void PhyloDistanceGamma< Tree>::setValue(const MatrixReal *dm, const MatrixReal *vm)
 {
	
	setDistanceMatrix(dm);
	setVarianceMatrix(vm);
 
 }
 */

void PhyloDistanceGamma::simulate( )
{
    this->redrawValue();
    
}

void PhyloDistanceGamma::setDistanceMatrix(const TypedDagNode< DistanceMatrix > *dm)
{
    
    // remove the old parameter first
    if ( distanceMatrix != NULL )
    {
        this->removeParameter( distanceMatrix );
        distanceMatrix = NULL;
    }
    
    // set the value
    distanceMatrix = dm;
    
    // add the new parameter
    this->addParameter( distanceMatrix );
    
    updateAlphaAndBetaMatrices( );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloDistanceGamma::setVarianceMatrix(const TypedDagNode< DistanceMatrix > *vm)
{
    
    // remove the old parameter first
    if ( varianceMatrix != NULL )
    {
        this->removeParameter( varianceMatrix );
        varianceMatrix = NULL;
    }
    
    // set the value
    varianceMatrix = vm;
    
    // add the new parameter
    this->addParameter( varianceMatrix );
    
    updateAlphaAndBetaMatrices( );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloDistanceGamma::setNames(const std::vector< std::string >& n)
{
    
    matrixNames = n;
    
}


/** Swap a parameter of the distribution */
void PhyloDistanceGamma::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
       // tau->getValue().getTreeChangeEventHandler().removeListener( this );
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
		//First, compute the distance matrix from the current tree
		
		//computeLnProbability();

		
      //  tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else if ( oldP == distanceMatrix )
    {
        
       // distanceMatrix = static_cast<const TypedDagNode<DistanceMatrix>* >( newP );
		setDistanceMatrix (static_cast<const TypedDagNode<DistanceMatrix>* >( newP ) );
		//computeLnProbability();
    }
    else if ( oldP == varianceMatrix )
    {
        
        //varianceMatrix = static_cast<const TypedDagNode<DistanceMatrix>* >( newP );
		setVarianceMatrix (static_cast<const TypedDagNode<DistanceMatrix>* >( newP ) );
		//computeLnProbability();
    }
    
}


void PhyloDistanceGamma::updateAlphaAndBetaMatrices( )
{
    //fill alpha and beta matrices
    // We assume the branch length is distributed according to a Gamma distribution.
    // We know the mean and the variance, we can compute the two parameters alpha and beta of the Gamma distribution.
    // mean = alpha/beta
    // variance = alpha /beta^2
    // Therefore:
    // alpha = mean^2 / variance
    // beta = mean / variance
    alphaMatrix = *( distanceMatrix->getValue().clone() );
    betaMatrix = *( varianceMatrix->getValue().clone() );
    double d, v;
    for (size_t i = 0; i < numTips - 1; ++i)
    {
        for (size_t j = i+1; j < numTips; ++j)
        {
            d = distanceMatrix->getValue()[i][j];
            v = varianceMatrix->getValue()[i][j];
            alphaMatrix[i][j] = alphaMatrix[j][i] =  pow(d, 2.0) / v;
            betaMatrix[i][j] = betaMatrix[j][i] = d / v ;
			std::cout <<"FILLING: alpha: "<< alphaMatrix[i][j] << " beta: " << betaMatrix[i][j] << "; mean: " << d << "; var: "<< v <<std::endl;
			
			
        }
    }
    
}
