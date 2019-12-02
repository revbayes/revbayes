// Example 4-3(c): An IsDerivedFrom constraints base 
// with testable value
//

#ifndef IsDerivedFrom_H
#define IsDerivedFrom_H

template<typename D, typename B>
class IsDerivedFrom
{
    class No { };
    class Yes { No no[3]; }; 
    
    static Yes Test( B* ); // not defined
    static No Test( ... ); // not defined 
    
    static void Constraints(D* p) { B* pb = p; pb = p; } 
    
public:
    enum { Is = sizeof(Test(static_cast<D*>(0))) == sizeof(Yes) }; 
    
    IsDerivedFrom() { void(*p)(D*) = Constraints; }
};

#endif
