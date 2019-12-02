//
//  PseudoGibbsMetaMove.h
//  
//
//  Created by Bastien Boussau on 05/07/13.
//  Copyright 2013 UC Berkeley. All rights reserved.
//
/**
 * @file
 * This file contains the declaration of PseudoGibbsMetaMove, which takes one or several moves,
 * and applies it/them independently to k subsets of the data. For each subset, the loglk difference 
 * between new and current state is computed. Among these k proposals, one is chosen randomly 
 * according to these differences.
 *
 * @brief Declaration of PseudoGibbsMetaMove
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-09-10 11:02:28 +0200 (Lun, 10 sep 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-09-08, version 1.0
 *
 * $Id: ScaleMove.h 1817 2012-09-10 09:02:28Z hoehna $
 */

#ifndef PseudoGibbsMetaMove_h
#define PseudoGibbsMetaMove_h


#include "Move.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include <vector>

namespace RevBayesCore {
    
    class PseudoGibbsMetaMove {
        
    public:
        ~PseudoGibbsMetaMove(void);                                             //!< Destructor
        
        // public method
        PseudoGibbsMetaMove*                            clone(void) const ;
        
    protected:
        PseudoGibbsMetaMove(const Move* m, AbstractHomologousDiscreteCharacterData* d, const size_t numOfSubsets);                                                                                                                                         //!< Default constructor
		PseudoGibbsMetaMove ( const PseudoGibbsMetaMove& p ) ;
		void buildSubsets();
		
        // Hidden member variables
        const Move*			                            move;
        size_t 											numberOfSubsets;
		AbstractHomologousDiscreteCharacterData*                  data;
		std::vector< std::vector<size_t> >				dataSubsets;
    };
    
}

#endif
