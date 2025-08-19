#ifndef UTIL_FILEFORMAT_H
#define UTIL_FILEFORMAT_H

#include "variant.h"

struct SeparatorFormat
{
    std::string separator;
    SeparatorFormat(const std::string& s):separator(s) {}
};

struct JSONFormat
{
};

/*
 * To add a new sample format 
 * 1. Create a struct FooFormat
 * 2. Add the struct to the SampleFormat std::variant type.
 */

typedef std::variant<SeparatorFormat,JSONFormat> SampleFormat;


#endif
