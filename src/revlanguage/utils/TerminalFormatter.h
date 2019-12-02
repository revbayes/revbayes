#ifndef TERMINAL_H
#define	TERMINAL_H

#include <string>

#ifndef TERMINAL_CODES
#define	TERMINAL_CODES
/**
 * Terminal codes: http://wiki.bash-hackers.org/scripting/terminalcodes
 */
namespace TerminalCodes {
    
#ifdef RB_XCODE
    const std::string _termReset("");
    const std::string _termUnderline("");
    const std::string _termBold("");
# else
    const std::string _termReset("\033[0m");
    const std::string _termUnderline("\033[4m");
    const std::string _termBold("\033[1m");
#endif
    
}

#endif

using namespace TerminalCodes;

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

