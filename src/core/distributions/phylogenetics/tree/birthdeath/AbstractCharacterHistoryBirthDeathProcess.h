#ifndef AbstractCharacterHistoryBirthDeathProcess_h
#define AbstractCharacterHistoryBirthDeathProcess_h


#include "Tree.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    class CharacterHistory;
    class DagNode;
    template <class valueType> class RbOrderedSet;
    
    class AbstractCharacterHistoryBirthDeathProcess : public TypedDistribution<Tree> {
        
    public:
        AbstractCharacterHistoryBirthDeathProcess(void);                                             //!< Constructor
        
        virtual                                            ~AbstractCharacterHistoryBirthDeathProcess(void);                //!< Virtual destructor
        
        // public member functions
        AbstractCharacterHistoryBirthDeathProcess*          clone(void) const = 0;                                          //!< Create an independent clone
        virtual double                                      computeLnProbability(void) = 0;                                 //!< Compute ln prob of current value
        virtual CharacterHistory&                           getCharacterHistory(void) = 0;                                      //!< Get the character histories
        virtual const CharacterHistory&                     getCharacterHistory(void) const = 0;                                      //!< Get the character histories
        virtual void                                        redrawValue(void) = 0;                                          //!< Draw a new random value from distribution
        virtual void                                        setValue(Tree *v, bool force) = 0;
        
    protected:
        // Parameter management functions
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter) = 0;                                      //!< get affected nodes
        virtual void                                        keepSpecialization(const DagNode* affecter) = 0;
        virtual void                                        restoreSpecialization(const DagNode *restorer) = 0;
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll) = 0;
        virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;            //!< Swap a
        
    private:
        
    };
    
}


#endif /* AbstractCharacterHistoryBirthDeathProcess_h */
