#include "MatrixRealSingleElementSwapProposal.h"

#include <cstddef>
#include <ostream>

#include "MatrixReal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;


MatrixRealSingleElementSwapProposal::MatrixRealSingleElementSwapProposal( StochasticNode<MatrixReal> *n, std::int64_t m, std::int64_t i ) : Proposal(),
    matrix( n ),
    margin( m ),
    index( i ),
    row_a( 0 ),
    col_a( 0 ),
    row_b( 0 ),
    col_b( 0 ),
    failed( false )
{
    addNode( matrix );
}


void MatrixRealSingleElementSwapProposal::cleanProposal( void )
{
}


MatrixRealSingleElementSwapProposal* MatrixRealSingleElementSwapProposal::clone( void ) const
{
    return new MatrixRealSingleElementSwapProposal( *this );
}


const std::string& MatrixRealSingleElementSwapProposal::getProposalName( void ) const
{
    static std::string name = "MatrixElementSwap";

    return name;
}


double MatrixRealSingleElementSwapProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Exchange two entries, anywhere (margin 0) or within one row (1) or column (2).
 */
double MatrixRealSingleElementSwapProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    MatrixReal& v = matrix->getValue();

    size_t n_rows = v.getNumberOfRows();
    size_t n_cols = v.getNumberOfColumns();

    // with no margin the whole matrix is one pool, otherwise the entries of a single line
    size_t n_lines   = ( margin == 0 ) ? 1 : ( ( margin == 1 ) ? n_rows : n_cols );
    size_t n_entries = ( margin == 0 ) ? n_rows * n_cols : ( ( margin == 1 ) ? n_cols : n_rows );

    if ( n_entries < 2 || n_lines == 0 )
    {
        failed = true;

        return 0.0;
    }

    failed = false;

    size_t line = 0;
    if ( margin != 0 )
    {
        if ( index < 0 )
        {
            line = size_t( rng->uniform01() * n_lines );
            if ( line >= n_lines ) line = n_lines - 1;
        }
        else
        {
            line = size_t( index );

            if ( line >= n_lines )
            {
                throw RbException() << "Cannot swap within line " << (line + 1) << " of a matrix with " << n_lines << " of them.";
            }
        }
    }

    size_t a = size_t( rng->uniform01() * n_entries );
    if ( a >= n_entries ) a = n_entries - 1;

    size_t b = size_t( rng->uniform01() * (n_entries - 1) );
    if ( b >= n_entries - 1 ) b = n_entries - 2;
    if ( b >= a ) b++;

    if ( margin == 0 )
    {
        row_a = a / n_cols;  col_a = a % n_cols;
        row_b = b / n_cols;  col_b = b % n_cols;
    }
    else if ( margin == 1 )
    {
        row_a = line;  col_a = a;
        row_b = line;  col_b = b;
    }
    else
    {
        row_a = a;  col_a = line;
        row_b = b;  col_b = line;
    }

    double tmp        = v[row_a][col_a];
    v[row_a][col_a]   = v[row_b][col_b];
    v[row_b][col_b]   = tmp;

    matrix->addTouchedElementIndex( row_a * n_cols + col_a );
    matrix->addTouchedElementIndex( row_b * n_cols + col_b );

    // a permutation of the value, so the proposal is symmetric
    return 0.0;
}


void MatrixRealSingleElementSwapProposal::prepareProposal( void )
{
}


void MatrixRealSingleElementSwapProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void MatrixRealSingleElementSwapProposal::undoProposal( void )
{
    if ( failed == false )
    {
        MatrixReal& v = matrix->getValue();

        double tmp      = v[row_a][col_a];
        v[row_a][col_a] = v[row_b][col_b];
        v[row_b][col_b] = tmp;
    }
}


void MatrixRealSingleElementSwapProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == matrix )
    {
        matrix = static_cast<StochasticNode<MatrixReal>* >(newN);
    }
}


void MatrixRealSingleElementSwapProposal::setProposalTuningParameter(double tp)
{
}


void MatrixRealSingleElementSwapProposal::tune( double rate )
{
}
