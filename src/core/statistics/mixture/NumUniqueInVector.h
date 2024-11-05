/**
 * @file
 * This file contains the declaration of the number of tables statistic for DPP-distributed variables, 
 * This simply computes the number of distinct value categories
 *
 * @brief Declaration of the NumUniqueInVector
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: NumUniqueInVector.h tracyh $
 */

#ifndef NumUniqueInVector_H
#define NumUniqueInVector_H

#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
	template <class valueType>
    class NumUniqueInVector : public TypedFunction<std::int64_t> {
        
    public:
																NumUniqueInVector(const TypedDagNode< RbVector<valueType> >* v);        //!< Default constructor
        virtual													~NumUniqueInVector(void);												//!< Destructor
        
        NumUniqueInVector&                                      operator=(const NumUniqueInVector& t);
        
        // Basic utility functions
        NumUniqueInVector*                                      clone(void) const;															//!< Clone object
        void													update(void);																//!< Clone the function
		
    protected:
        void													swapParameterInternal(const DagNode *oldP, const DagNode *newP);			//!< Implementation of swaping parameters
        size_t												    findValueinValuePerTable(valueType v, std::vector<valueType> vpt);
    private:
		
		// members
        const TypedDagNode< RbVector<valueType> >*              elementVals;
    };
	
}

using namespace RevBayesCore;

template<class valueType>
NumUniqueInVector<valueType>::NumUniqueInVector(const TypedDagNode< RbVector<valueType> >* v) : TypedFunction<std::int64_t>( new long(0) ), elementVals( v ) {
    // add the tree parameter as a parent
    addParameter( elementVals );
    update();
}

template<class valueType>
NumUniqueInVector<valueType>::~NumUniqueInVector( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

template<class valueType>
NumUniqueInVector<valueType>* NumUniqueInVector<valueType>::clone( void ) const {
    return new NumUniqueInVector<valueType>( *this );
}

template<class valueType>
void NumUniqueInVector<valueType>::update( void ) {
    
	const RbVector<valueType>& pv = elementVals->getValue();
	std::vector<valueType> valuePerTable;
	std::vector<int> numCustomerPerTable;
	int numTables = 0;
	size_t numElements = pv.size();
	for (int i=0; i<numElements; i++){
		valueType v = pv[i];
		size_t tID = findValueinValuePerTable(v, valuePerTable);
		if (tID < valuePerTable.size()){
			numCustomerPerTable[tID] += 1;
		}
		else{
			valuePerTable.push_back(v);
			numCustomerPerTable.push_back(1);
			numTables++;
		}
	}


	*value = numTables;
}

template<class valueType>
void NumUniqueInVector<valueType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    if (oldP == elementVals) {
        elementVals = static_cast< const TypedDagNode< RbVector<valueType> >* >( newP );
    }
}

template <class valueType>
size_t NumUniqueInVector<valueType>::findValueinValuePerTable(valueType v, std::vector<valueType> vpt){
	
	size_t tID = vpt.size();
	if (tID == 0)
		return tID;
	else{
		for (size_t i=0; i<vpt.size(); i++){
			if (vpt[i] == v){
				return i;
			}
		}
	}
	return tID;
}


#endif

