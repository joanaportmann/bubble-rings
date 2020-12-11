// Use this to enable floating point exceptions only in a single scope (RAII)
//
//
// float x = 1./0.;
// {
//     FloatExceptionEnabler foo;
//     x = 2./0.;
// }
// // here float exceptions are restored to the 'before' state

#pragma once

#include <cfenv> // NOLINT(build/c++11)
#include "assert.hh"

class FloatExceptionEnabler {
public:
    explicit FloatExceptionEnabler(std::fexcept_t enable = FE_DIVBYZERO
                                                | FE_INVALID
                                                | FE_OVERFLOW
                                                | FE_UNDERFLOW)
    {
        release_assert(0 == std::fegetexceptflag(&saved_flags_, FE_ALL_EXCEPT));
        release_assert(0 == std::fesetexceptflag(&enable, FE_ALL_EXCEPT));
    }
    ~FloatExceptionEnabler() {
        release_assert(0 == std::fesetexceptflag(&saved_flags_, FE_ALL_EXCEPT));
    }
private:
    std::fexcept_t saved_flags_;
};
