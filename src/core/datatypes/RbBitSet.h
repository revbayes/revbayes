#ifndef RbRbBitSet_H
#define RbRbBitSet_H

#include <stddef.h>
#include <ostream>
#include <vector>

namespace RevBayesCore {
    
    /**
     * RevBayes class for bit sets.
     *
     * This class is essentially a wrapper of a std::vector<bool> with some additional functional.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-09-06, version 1.0
     */
    class RbBitSet {
        
    public:
        RbBitSet(void);                                                                                         //!< Constructor requires character type
        RbBitSet(size_t n, bool def = false);                                                                   //!< Constructor requires character type
        virtual                        ~RbBitSet(void) {}
        
        bool                            operator[](size_t i) const;

        bool                            operator==(const RbBitSet &bs) const;
        bool                            operator!=(const RbBitSet &bs) const;
        bool                            operator<(const RbBitSet &bs) const;
        
        RbBitSet                        operator&(const RbBitSet &bs) const;
        RbBitSet                        operator|(const RbBitSet &bs) const;
        RbBitSet                        operator^(const RbBitSet &bs) const;
        RbBitSet&                       operator~();
        RbBitSet&                       operator&=(const RbBitSet &bs);
        RbBitSet&                       operator|=(const RbBitSet &bs);


        void                            clear(void);
        bool                            empty(void) const;
        void                            flip();
        void                            flip(size_t i);
        size_t                          getNumberSetBits(void) const;                                           //!< Get the number of bits set.
        size_t                          getFirstSetBit(void) const;                                             //!< Get the number of bits set.
        bool                            isSet(size_t i) const;
        void                            resize(size_t size);
        void                            set(size_t i);
        size_t                          size(void) const;
        void                            unset(size_t i);

        
    private:

        std::vector<bool>               value;
        size_t                          num_set_bits;

        
    };
    
    // Global functions using the class
    std::ostream&                               operator<<(std::ostream& o, const RbBitSet& x);                    //!< Overloaded output operator
    
}

#endif
