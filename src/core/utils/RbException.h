#ifndef RbException_H
#define RbException_H


#include <iostream>


class RbException {

    public:
        // Exception types
        enum                        ExceptionType { DEFAULT,
                                                    MATH_ERROR,
                                                    MISSING_VARIABLE,
                                                    QUIT };         //!< Exception types
        static std::string          ExceptionName[];                                        //!< Exception type names

                                    RbException(void);                                      //!< Default constructor
                                    RbException(const std::string& msg);                    //!< Default with message 
                                    RbException(ExceptionType type, const std::string& msg="");//!< General constructor

        // Regular functions
        ExceptionType               getExceptionType(void) const;                           //!< Get exception type
        void                        setMessage(std::string msg);
        std::string                 getMessage(void) const;
        void                        print(std::ostream &o) const;                           //!< Print the exception
     
    private:
	    ExceptionType               exception_type;                                         //!< Exception type
	    std::string                 message;                                                //!< Error message
    
};

#endif

