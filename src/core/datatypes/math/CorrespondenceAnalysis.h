//
// File: CorrespondenceAnalysis.h
// Created by: Mathieu Groussin
// Created on: Sun Mar 06 10:03 2011
//

/*
 Copyright or ¬© or Copr. Bio++ Development Team, (November 17, 2004)
 
 This software is a computer program whose purpose is to provide basal and
 utilitary classes. This file belongs to the Bio++ Project.
 
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


#ifndef _CORRESPONDENCEANALYSIS_H_
#define _CORRESPONDENCEANALYSIS_H_

#include "DualityDiagram.h"
#include "MatrixReal.h"
#include "RbException.h"

namespace RevBayesCore
/**
 * @brief This class allows to perform a correspondence analysis.
 *
 * All values in the input table have to be non-negative.
 * The DualityDiagram class, core class of a multivariate analysis, is called internally.
 *
 * The code of this class is deeply inspired from the R code of the dudi.coa function available in the ade4 package.
 */
{
    class CorrespondenceAnalysis:
    public DualityDiagram
    {
        
    public:
        /**
         * @brief Build a new CorrespondenceAnalysis object.
         *
         * @param data The input data (a RowMatrix) to analyse.
         * @param nbAxes The number of kept axes during the analysis.
         * @param tol Tolerance threshold for null eigenvalues (a value less than tol times the first one is considered as null)
         * @param verbose Should warnings be dispayed.
         * @throw Exception if an error occured.
         */
        CorrespondenceAnalysis(
                               const MatrixReal& data,
                               unsigned int nbAxes,
                               double tol);
        
        virtual ~CorrespondenceAnalysis() {}
        
        CorrespondenceAnalysis* clone() const { return new CorrespondenceAnalysis(*this); }
        
    public:
        double getSumOfAllValues() const { return n_; }
    
    private:
        double n_;

    };
} // end of namespace RevBayesCore.

#endif  // _CORRESPONDENCEANALYSIS_H_

