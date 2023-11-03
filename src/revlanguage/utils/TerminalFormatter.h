#ifndef TERMINAL_H
#define	TERMINAL_H

#include <string>

#ifndef TERMINAL_CODES
#define	TERMINAL_CODES
/**
 * Terminal codes: http://wiki.bash-hackers.org/scripting/terminalcodes
 */
namespace TerminalCodes
{
    extern std::string _termReset;
    extern std::string _termUnderline;
    extern std::string _termBold;
}

#endif

using namespace TerminalCodes;

void enableTermAnsi();

class TerminalFormatter {
public:
    
    static std::string makeBold(std::string s) {
        return _termBold + s + _termReset;
    }
    
    static std::string makeUnderlined(std::string s) {
        return _termUnderline + s + _termReset;
    }
};


#endif	/* TERMINAL_H */

