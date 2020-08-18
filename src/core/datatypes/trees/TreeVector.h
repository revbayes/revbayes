#ifndef TreeVector_h
#define TreeVector_h

#include <stdio.h>
#include <vector>
#include <iosfwd>

#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {


    /**
     * @brief Vector of trees
     *
     *
     * @copyright Copyright 2009-
     * @author
     * @since 2014-03-18, version 1.0
     */
    class TreeVector : public RbVector<Tree>, public MemberObject<Tree> {

    public:

        TreeVector(void);                                              //!< Default constructor
        virtual                             ~TreeVector(void);

        bool                                operator==(const TreeVector &mve) const;
        bool                                operator!=(const TreeVector &mve) const;

        // public methods
        void                                addValues(const Tree *t, long n);
        void                                clear(void);
        TreeVector*                         clone(void) const;
        size_t                              getNumberOfValues(void) const;
        std::vector<Tree>&                  getValues(void);                                                //!< Get the values for this element.
        const std::vector<Tree>&            getValues(void) const;                                          //!< Get the values for this element.
        void                                setNumberOfTrees(long n);
        void                                setValues(const std::vector<Tree> &v, const long n);

        // Utility functions
        void                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Tree &rv) const;


    private:

        // private members
        long                                num_trees;
        std::vector<Tree>                   values;

    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const TreeVector& x);                                         //!< Overloaded output operator
}


#endif /* TreeVector_hpp */
