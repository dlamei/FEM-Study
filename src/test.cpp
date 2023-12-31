#include "utils.h"
#include <thread>
#include <csignal>

void passed(const char *test_name) {
    std::cout << "[";
    std::cout << "\033[32m" << "PASSED" << "\033[0m";
    std::cout << "]\t" << test_name << "\n";
}

void failed(const char *test_name, const char *msg) {
    std::cout << "[";
    std::cout << "\033[1;31m" << "FAILED" << "\033[0m";
    std::cout << "]\t" << test_name << "\n";
    std::cout << "width message:\n" << msg << "\n";
}

void run_all_tests() {
    auto &tests = global_test_collector::get_tests();

    for (auto test : tests) {
        auto msg = test.fn();
        if (msg) {
            failed(test.name, msg->c_str());
        } else {
            passed(test.name);
        }
    }
}

#ifdef COMPILE_TESTS

int main() {
    run_all_tests();
}

#endif
