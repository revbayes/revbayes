//
//  ComputeLtFunction.cpp
//
//  Created by Rachel Warnock, Marc Manceau 30.01.2020.
//
#include "ComputeLtFunction.h"

#include <vector>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "TypedDagNode.h"
// other includes

namespace RevBayesCore { class DagNode; } //tt

using namespace RevBayesCore;

/** ComputeLtFunction
 * this is temporary until we work out where this function should live in the long term
 * @param v the vector of branching times
 */
ComputeLtFunction::ComputeLtFunction(
	const TypedDagNode< RbVector<double> > *a,
	const TypedDagNode< RbVector<double> > *b,
	const TypedDagNode< RbVector<double> > *c,
	const TypedDagNode< RbVector<double> > *d,
	const TypedDagNode< RbVector<double> > *e,
	const TypedDagNode< RbVector<double> > *f,
	const TypedDagNode< double > *t,
	const TypedDagNode< double > *l,
	const TypedDagNode< double > *m,
	const TypedDagNode< double > *p,
	const TypedDagNode< double > *o,
	const TypedDagNode< double > *rho,
	const TypedDagNode< double > *r
 	) : TypedFunction<double>( new double(0.0) ),
	listA ( a ), listB ( b ), listC ( c ), listD ( d ), listE ( e ), listF ( f ),
	tor ( t ), lambda ( l ), mu ( m ), psi ( p ), omega ( o ), rho ( rho ), removalPr ( r )
{

	//add the parameter as parents
	this->addParameter( listA );
	this->addParameter( listB );
	this->addParameter( listC );
	this->addParameter( listD );
	this->addParameter( listE );
	this->addParameter( listF );

	this->addParameter( tor );
	this->addParameter( lambda );
	this->addParameter( mu );
	this->addParameter( psi );
	this->addParameter( omega );
	this->addParameter( rho );
	this->addParameter( removalPr ); // I don't think you need to do this

	update();

}

ComputeLtFunction::~ComputeLtFunction( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/** clone function */
ComputeLtFunction* ComputeLtFunction::clone( void ) const
{
    return new ComputeLtFunction( *this );
}

// required -> update function, where you do the work
void ComputeLtFunction::update( void ) {

	const std::vector<double> &a = listA->getValue(); // branching times
	const std::vector<double> &b = listB->getValue(); // sampled ancestors
	const std::vector<double> &c = listC->getValue(); // terminal non-removed
	const std::vector<double> &d = listD->getValue(); // terminal removed
	const std::vector<double> &e = listE->getValue(); //
	const std::vector<double> &f = listF->getValue();

	double out = 0.0;
	for ( std::vector<double>::const_iterator it = a.begin(); it != a.end(); ++it) {
			out += *it;
	}

	out = lambda->getValue() * out;

	MatrixReal B(2, 2, out);

	RbMath::expMatrixPade(B, B, 4);

	*this->value = B[0][0]; // this will eventually be the lk of the tree + occurrences

}

// required -> swap function. See Variance Function for inspiration
void ComputeLtFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP ) {

	if ( oldP == listA	)
	{
			listA = static_cast<const TypedDagNode< RbVector<double> >* >( newP ); //tt
	}

}
