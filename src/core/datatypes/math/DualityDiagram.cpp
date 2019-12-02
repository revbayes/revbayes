//
// File: DualityDiagram.cpp
// Created by: Mathieu Groussin
//

/*
 Copyright or ¬© or Copr. Bio++ Development Tools, (November 16, 2004)
 
 This software is a computer program whose purpose is to provide classes
 for phylogenetic data analysis.
 
 This software is governed by the CeCILL  license under French law and
 abiding by the rules of distribution of free software.  You can  use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".
 
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability.
 
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and,  more generally, to use and operate it in the
 same conditions as regards security.
 
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.
 */

#include "DualityDiagram.h"

#include <cmath>

#include "EigenSystem.h"
#include "MatrixReal.h"
#include "RbMathMatrix.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

DualityDiagram::DualityDiagram(
                               const MatrixReal& matrix,
                               const std::vector<double>& rowWeights,
                               const std::vector<double>& colWeights,
                               unsigned int nbAxes,
                               double tol) :
rowWeights_(rowWeights),
colWeights_(colWeights),
nbAxes_(nbAxes),
eigenValues_(),
eigenVectors_(2,2),
rowCoord_(2,2),
colCoord_(2,2),
ppalAxes_(2,2),
ppalComponents_(2,2)
{
    check_(matrix, rowWeights, colWeights, nbAxes);
    compute_(matrix, tol);
}

void DualityDiagram::check_(
                            const MatrixReal& matrix,
                            const std::vector<double>& rowWeights,
                            const std::vector<double>& colWeights,
                            unsigned int nbAxes)
{
    size_t rowNb = matrix.getNumberOfRows();
    size_t colNb = matrix.getNumberOfColumns();
    
    if (rowWeights.size() != rowNb)
        throw RbException("DualityDiagram::check_. The number of row weigths has to be equal to the number of rows!");
    if (colWeights.size() != colNb)
        throw RbException("DualityDiagram::check_. The number of column weigths has to be equal to the number of columns!");
    
    // All row weigths have to be positive
    for (std::vector<double>::const_iterator it = rowWeights.begin(); it != rowWeights.end(); it++)
    {
        if (*it < 0.)
            throw RbException("DualityDiagram::check_. All row weights have to be positive");
    }
    
    // All column weigths have to be positive
    for (std::vector<double>::const_iterator it = colWeights.begin(); it != colWeights.end(); it++)
    {
        if (*it < 0.)
            throw RbException("DualityDiagram::check_. All column weights have to be positive");
    }
}


void DualityDiagram::setData(
                             const MatrixReal& matrix,
                             const std::vector<double>& rowWeights,
                             const std::vector<double>& colWeights,
                             unsigned int nbAxes,
                             double tol)
{
    check_(matrix, rowWeights, colWeights, nbAxes);
    rowWeights_ = rowWeights;
    colWeights_ = colWeights;
    nbAxes_     = nbAxes;
    compute_(matrix, tol);
}

void DualityDiagram::compute_(const MatrixReal& matrix, double tol)
{
    size_t rowNb = matrix.getNumberOfRows();
    size_t colNb = matrix.getNumberOfColumns();
    
    // If there are less rows than columns, the variance-covariance or correlation matrix is obtain differently (see below)
    bool transpose = (rowNb < colNb);
    
    // The initial matrix is multiplied by the square root of the row weigths.
    std::vector<double> rW(rowWeights_);
    for (size_t i = 0; i < rowWeights_.size(); i++)
    {
        rW[i] = sqrt(rowWeights_[i]);
    }
    
    MatrixReal M1(rowNb, colNb);
    RbMath::hadamardMult(matrix, rW, M1, true);
    
    // The resulting matrix is then multiplied by the square root of the column weigths.
    std::vector<double> cW(colWeights_);
    for (unsigned int i = 0; i < colWeights_.size(); i++)
    {
        cW[i] = sqrt(colWeights_[i]);
    }
    
    MatrixReal M2(rowNb,colNb);
    RbMath::hadamardMult(M1, cW, M2, false);
    
    // The variance-covariance (if the data is centered) or the correlation (if the data is centered and normalized) matrix is calculated
    MatrixReal tM2(M2.getNumberOfColumns(),M2.getNumberOfRows());
    RbMath::transposeMatrix(M2, tM2);
    MatrixReal *M3 = new MatrixReal(0,0);
    if (!transpose)
        (*M3) = tM2 * M2;
    else
        (*M3) = M2 * tM2;
    
    EigenSystem eigen(M3);
    eigen.update();
    // @todo: This may be implemented some time ... (Sebastian)
//    if (!eigen.isSymmetric())
//        throw RbException("DualityDiagram (constructor). The variance-covariance or correlation matrix should be symmetric...");
    
    eigenValues_ = eigen.getRealEigenvalues();
    eigenVectors_ = eigen.getEigenvectors();
    
    // How many significant axes have to be conserved?
    size_t rank = 0;
    for (size_t i = 0; i < eigenValues_.size(); ++i)
    {
        if ((eigenValues_[i] / eigenValues_[0]) > tol)
            rank++;
    }
    
    if (nbAxes_ <=0)
    {
        throw RbException("DualityDiagram (constructor). The number of axes to keep must be positive.");
    }
    if (nbAxes_ > rank)
    {
        nbAxes_ = rank;
    }
    
//    /*The eigen values are initially sorted into ascending order by the 'eigen' function. Here the significant values are sorted
//     in the other way around.*/
//    std::vector<double> tmpEigenValues(nbAxes_);
//    size_t cpt = 0;
//    for (size_t i = eigenValues_.size(); i > (eigenValues_.size() - nbAxes_); i--)
//    {
//        tmpEigenValues[cpt] = eigenValues_[i-1];
//        cpt++;
//    }
//    eigenValues_ = tmpEigenValues;
    
    for (std::vector<double>::iterator it = rowWeights_.begin(); it != rowWeights_.end(); it++)
    {
        if (*it == 0.)
            *it = 1.;
    }
    
    for (std::vector<double>::iterator it = colWeights_.begin(); it != colWeights_.end(); it++)
    {
        if (*it == 0.)
            *it = 1.;
    }
    
    std::vector<double> dval(nbAxes_);
    for (size_t i = 0; i < dval.size(); i++)
    {
        dval[i] = sqrt(eigenValues_[i]);
    }
    
    std::vector<double> invDval(nbAxes_);
    for (size_t i = 0; i < invDval.size(); i++)
    {
        invDval[i] = 1. / sqrt(eigenValues_[i]);
    }
    
    // Calculation of the row and column coordinates as well as the principal axes and components:
    if (!transpose)
    {
        std::vector<double> tmpColWeights(colNb);
        for (unsigned int i = 0; i < colWeights_.size(); i++)
        {
            tmpColWeights[i] = 1. / sqrt(colWeights_[i]);
        }
        
        // The eigen vectors are placed in the same order as their corresponding eigen value in eigenValues_.
        MatrixReal tmpEigenVectors(0,0);
        tmpEigenVectors.resize(eigenVectors_.getNumberOfRows(), nbAxes_);
        size_t cpt2 = 0;
        for (size_t i = 0; i < nbAxes_; i++)
        {
            for (unsigned int j = 0; j < eigenVectors_.getNumberOfRows(); j++)
            {
                tmpEigenVectors[j][cpt2] = eigenVectors_[j][i];
            }
            cpt2++;
        }

        
        // matrix of principal axes
        RbMath::hadamardMult(tmpEigenVectors, tmpColWeights, ppalAxes_, true);
        // matrix of row coordinates
        MatrixReal tmpRowCoord_(0,0);
        tmpRowCoord_.resize(rowNb, nbAxes_);
        RbMath::hadamardMult(matrix, colWeights_, tmpRowCoord_, false);
        rowCoord_ = tmpRowCoord_ * ppalAxes_;
        
        // matrix of column coordinates
        RbMath::hadamardMult(ppalAxes_, dval, colCoord_, false);
        // matrix of principal components
        RbMath::hadamardMult(rowCoord_, invDval, ppalComponents_, false);
    }
    else
    {
        std::vector<double> tmpRowWeights(rowNb);
        for (unsigned int i = 0; i < rowWeights_.size(); i++)
        {
            tmpRowWeights[i] = 1. / sqrt(rowWeights_[i]);
        }
        
        // The eigen vectors are placed in the same order as their corresponding eigen value in eigenValues_.
        MatrixReal tmpEigenVectors(0,0);
        tmpEigenVectors.resize(eigenVectors_.getNumberOfRows(), nbAxes_);
        size_t cpt2 = 0;
        for (size_t i = 0; i < nbAxes_; i++)
        {
            for (unsigned int j = 0; j < eigenVectors_.getNumberOfRows(); j++)
            {
                tmpEigenVectors[j][cpt2] = eigenVectors_[j][i];
            }
            cpt2++;
        }

        
        // matrix of principal components
        RbMath::hadamardMult(tmpEigenVectors, tmpRowWeights, ppalComponents_, true);
        // matrix of column coordinates
        MatrixReal tmpColCoord_(colNb, nbAxes_);
        RbMath::hadamardMult(matrix, rowWeights_, tmpColCoord_, true);
        MatrixReal tTmpColCoord_(tmpColCoord_.getNumberOfColumns(),tmpColCoord_.getNumberOfRows());
        RbMath::transposeMatrix(tmpColCoord_, tTmpColCoord_);
        colCoord_ = tTmpColCoord_ * ppalComponents_;
        
        // matrix of row coordinates
        RbMath::hadamardMult(ppalComponents_, dval, rowCoord_, false);
        // matrix of principal axes
        RbMath::hadamardMult(colCoord_, invDval, ppalAxes_, false);
    }
}

/******************************************************************************/

DualityDiagram::~DualityDiagram() {}

/******************************************************************************/

