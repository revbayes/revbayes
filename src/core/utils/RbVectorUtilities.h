/**
 * @file
 * This file helper functions for manipulating RlStrings in RevBayes.
 *
 * @brief Namespace containing helper functions for RlStrings
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */
#ifndef RlVectorUtilities_H
#define RlVectorUtilities_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace VectorUtilities {
       
    /**
     * Generic to_string function
     * @param value
     * @return
     */
    template <typename T>
    std::string to_string(T value) {
        std::ostringstream os;
        os << value;
        return os.str();
    }
    
    /**
     * Find and return the index of pivot element.
     * @param a - The array.
     * @param first - The start of the sequence.
     * @param last - The end of the sequence.
     * @return - the pivot element
     */
    template <typename T>
    int pivot(std::vector<T> &v, int first, int last)
    {
        int  p = first;
        const T& pivotElement = v[first];
        
        for (int i = first+1 ; i <= last ; i++)
        {
            /* If you want to sort the list in the other order, change "<=" to ">" */
            if ( v[i] <= pivotElement)
            {
                p++;
                
                T temp = T( v[i] );
                v[i] = v[p];
                v[p] = temp;
            }
        }
        
        
        T temp = T( v[p] );
        v[p] = v[first];
        v[first] = temp;
        
        return p;
    }
    
    
    template <class T>
    void quicksort(std::vector<T> &v, int first, int last)
    {
        if (first < last)
        {
            int pivotElement = pivot<T>(v, first, last);
            quicksort<T>(v, first, pivotElement-1);
            quicksort<T>(v, pivotElement+1, last);
        }
    }
    
    
    template <class T>
    void sort(std::vector<T> &v)
    {
        // just delegate to our internal quicksort method.
        quicksort<T>(v, 0, int(v.size())-1);
    }

    template <typename T>
    T sum(const std::vector<T> & v)
    {
	T total = 0;
	for(int i=0; i < v.size(); i++)
	    total += v[i];
	return total;
    }

    template <typename T>
    void scale(std::vector<T> & v, T factor)
    {
	for(int i=0; i < v.size(); i++)
	    v[i] *= factor;
    }

    inline void normalize(std::vector<double> & v)
    {
	scale(v, 1.0/sum(v));
    }
}


#endif
