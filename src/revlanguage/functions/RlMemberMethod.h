#ifndef RlMemberMethod_H
#define RlMemberMethod_H


namespace RevLanguage {
    
    
    class Argument;
    
    /**
     * @brief Rev interface for member methods.
     *
     * This interface specifies that classes implementing the interface will be a member methods,
     * either function or procedure.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-10-29
     *
     */
    class MemberMethod {
        
    public:
        
        virtual void            setMemberObject(const RevPtr<RevVariable> &obj) = 0;           //!< Set the member object to which this method belongs
	virtual ~MemberMethod() {};
    };
    
}


#endif

