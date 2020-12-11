
#include "assert.hh"

#include <array>
#include <cstdlib>
#include <iostream>
#include <mutex> // NOLINT(build/c++11)

#ifndef _MSC_VER
extern "C" {
#include <execinfo.h>
}
#endif

namespace {
std::mutex assert_mutex;
} // namespace
void __release_assert_fail(const char *__assertion, const char *__file,
                           unsigned int __line, const char *__function,
                           const std::string &extra)
{
    std::lock_guard<std::mutex> lock(assert_mutex);
    std::cout << std::flush;
    std::cerr << std::flush;

    std::cerr << "\n\nAssertion failed (" << __assertion
              << ") in function " << __function
              << " in file " << __file
              << ":" << __line;
    if (!extra.empty()) {
        std::cerr << " (" << extra << ")";
    }
    std::cerr << std::endl;
#ifndef _MSC_VER
    std::array<void*, 100> bt_buf;
    int size = backtrace(&bt_buf.front(), static_cast<int>(bt_buf.size()));
    backtrace_symbols_fd(&bt_buf.front(), size, 2);
#endif
    std::abort();
}
