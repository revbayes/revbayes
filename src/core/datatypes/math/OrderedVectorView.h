#ifndef OrderedVectorView_H
#define OrderedVectorView_H

#include "OrderedVector.h"

#include <cstddef>

namespace RevBayesCore {

    /**
     * @brief A value seen as the ordered vectors it contains.
     *
     * Not a container: the entries stay where the value keeps them, and this presents them a
     * chain at a time so a move can find an entry's neighbours. Whether a chain is a row, a
     * column or something with a stride is the value's business, so a move never branches on
     * the shape of what it was handed.
     *
     * The chains are totally ordered, deliberately. Moves written against OrderedVector assume
     * each entry has one neighbour on either side, which an arbitrary partial order would not give.
     */
    class OrderedVectorView {

    public:
        virtual ~OrderedVectorView(void) {}

        virtual size_t          numOrderedVectors(void) const = 0;
        virtual OrderedVector   getOrderedVector(size_t k) = 0;
    };

}

#endif
