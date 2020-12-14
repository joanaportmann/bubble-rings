#pragma once
//
// Created by Martin Heistermann on 26.07.18.
//

#include <string>
#include <cassert>


#ifdef _MSC_VER
#  define NORETURN(x) __declspec(noreturn) x
#  define FUNCNAME __FUNCSIG__
#else
#  define NORETURN(x) x __attribute__((__noreturn__))
#  define FUNCNAME __PRETTY_FUNCTION__
#endif

# define release_assert_extra(expr, extra) \
  ((expr) \
   ? static_cast<void> (0)	\
   : __release_assert_fail (#expr, __FILE__, __LINE__, FUNCNAME, extra))

# define release_assert(expr) release_assert_extra(expr, "")

#ifndef NDEBUG
# define debug_assert_extra(expr, extra) release_assert_extra(expr, extra)
#else
# define debug_assert_extra(expr, extra) static_cast<void> (0)
#endif

# define debug_assert(expr) debug_assert_extra(expr, "")


// #ifdef assert
// #undef assert
// #endif

//#define assert(x) debug_assert(x)

NORETURN(void __release_assert_fail (const char *__assertion, const char *__file,
                 unsigned int __line, const char *__function, const std::string &extra));

