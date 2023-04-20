#include "TerminalFormatter.h"

namespace TerminalCodes
{
    std::string _termReset;
    std::string _termUnderline;
    std::string _termBold;
}

void enableTermAnsi()
{
    _termReset = "\033[0m";
    _termUnderline = "\033[4m";
    _termBold = "\033[1m";
}

