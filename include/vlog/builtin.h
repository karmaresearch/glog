#ifndef _BUILTIN_H
#define _BUILTIN_H

#include <vlog/term.h>

#include <functional>
#include <array>

struct BuiltinFunction {
    //Max five args
    std::array<uint8_t, 5> posArgs;
    std::function<bool(Term_t * , uint8_t *)> fn;
};

#endif
