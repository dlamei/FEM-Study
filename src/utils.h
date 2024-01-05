#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <sstream>
#include <optional>

#include <cstddef>
#include <stdint.h>
typedef uint8_t   u8;
typedef  int8_t   i8;
typedef uint16_t  u16;
typedef  int16_t  i16;
typedef uint32_t  u32;
typedef  int32_t  i32;
typedef uint64_t  u64;
typedef  int64_t  i64;
typedef size_t    usize;
typedef ptrdiff_t isize;
typedef float     f32;
typedef double    f64;

static_assert(sizeof(u8) == sizeof(i8));
static_assert(sizeof(u16) == sizeof(i16));
static_assert(sizeof(u32) == sizeof(i32));
static_assert(sizeof(u64) == sizeof(i64));
static_assert(sizeof(usize) == sizeof(isize));

static_assert(sizeof(u8) == 1);
static_assert(sizeof(u16) == 2);
static_assert(sizeof(u32) == 4);
static_assert(sizeof(f32) == 4);
static_assert(sizeof(u64) == 8);
static_assert(sizeof(f64) == 8);

typedef f32 scalar;

#define SCALAR_EPS std::numeric_limits<scalar>::epsilon()

inline bool cmp_scalar(scalar a, scalar b, scalar eps = SCALAR_EPS) {
	return std::fabs(a - b) <= eps;
}

template<typename T, typename = void>
constexpr bool is_defined = false;

template<typename T>
constexpr bool is_defined<T, decltype(typeid(T), void())> = true;


/* defer for c++ ( from [https://github.com/gingerBill/gb] )*/

template <typename T> struct gbRemoveReference       { typedef T Type; };
template <typename T> struct gbRemoveReference<T &>  { typedef T Type; };
template <typename T> struct gbRemoveReference<T &&> { typedef T Type; };

template <typename T> inline T &&gb_forward(typename gbRemoveReference<T>::Type &t)  { return static_cast<T &&>(t); }
template <typename T> inline T &&gb_forward(typename gbRemoveReference<T>::Type &&t) { return static_cast<T &&>(t); }
template <typename T> inline T &&gb_move   (T &&t)                                   { return static_cast<typename gbRemoveReference<T>::Type &&>(t); }
template <typename F>
struct gbprivDefer {
    F f;
    gbprivDefer(F &&f) : f(gb_forward<F>(f)) {}
    ~gbprivDefer() { f(); }
};
template <typename F> gbprivDefer<F> gb__defer_func(F &&f) { return gbprivDefer<F>(gb_forward<F>(f)); }


/* assert with error messages */

inline void gb__assert(const char* expr_str, bool expr, const char* file, const char *func, int line, const char* msg = nullptr)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t";
        if (msg) std::cerr << msg;
        std::cerr << "\nFailed:\t\t" << expr_str << "\n";
        std::cerr << "Function:\t" << func << "\n";
        std::cerr << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

#define GB_DEFER_1(x, y) x##y
#define GB_DEFER_2(x, y) GB_DEFER_1(x, y)
#define GB_DEFER_3(x)    GB_DEFER_2(x, __COUNTER__)
#define defer(code)      auto GB_DEFER_3(_defer_) = gb__defer_func([&]()->void{code;})

// from [boost/current_function.hpp](https://www.boost.org/doc/libs/1_62_0/boost/current_function.hpp)
#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600)) || defined(__ghs__)

# define BOOST_CURRENT_FUNCTION __PRETTY_FUNCTION__

#elif defined(__DMC__) && (__DMC__ >= 0x810)

# define BOOST_CURRENT_FUNCTION __PRETTY_FUNCTION__

#elif defined(__FUNCSIG__)

# define BOOST_CURRENT_FUNCTION __FUNCSIG__

#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))

# define BOOST_CURRENT_FUNCTION __FUNCTION__

#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)

# define BOOST_CURRENT_FUNCTION __FUNC__

#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)

# define BOOST_CURRENT_FUNCTION __func__

#elif defined(__cplusplus) && (__cplusplus >= 201103)

# define BOOST_CURRENT_FUNCTION __func__

#else

# define BOOST_CURRENT_FUNCTION "(unknown)"

#endif

// because macros are too hard for microsoft
// https://developercommunity.visualstudio.com/t/-va-args-does-not-work-with-multiple-macros/167330k
#define EXPAND_VA_ARGS( x ) x

#define GB_ASSERT_1(x) gb__assert(#x, x, __FILE__, BOOST_CURRENT_FUNCTION, __LINE__, nullptr)
#define GB_ASSERT_2(x, msg) gb__assert(#x, x, __FILE__, BOOST_CURRENT_FUNCTION, __LINE__, msg)
#define GET_GB_ASSERT(_1, _2, NAME, ...) NAME
#ifdef assert
#undef assert
#endif
#define assert(...) EXPAND_VA_ARGS( GET_GB_ASSERT(__VA_ARGS__, GB_ASSERT_2, GB_ASSERT_1) (__VA_ARGS__) )

#ifdef NDEBUG
#define db_assert(...)
#else
#define db_assert(...) assert(__VA_ARGS__)
#endif

typedef std::optional<std::string>(* test_fn)();

struct global_test_collector {

    struct Test {
        test_fn fn;
        const char *name;
    };

    static global_test_collector &get_inst() {
        static global_test_collector instance{};
        return instance;
    }

    static const std::vector<Test> &get_tests() {
        return global_test_collector::get_inst().tests;
    }

    static void push(test_fn test, const char *name) {
        global_test_collector::get_inst().tests.push_back(Test {.fn = test, .name = name });
    }

    global_test_collector(global_test_collector const &) = delete;
    void operator=(global_test_collector const &) = delete;

private:
    global_test_collector() = default;
    std::vector<Test> tests{};
};


struct __register_test {
    explicit __register_test(test_fn test, const char *test_name) {
        global_test_collector::push(test, test_name);
    }
};

#ifdef COMPILE_TESTS

#define TEST(test_name, test_body) \
    std::optional<std::string> test_##test_name() test_body \
    __register_test register_##test_name(test_##test_name, #test_name);


inline std::optional<std::string> __test_assert(const char* expr_str, bool expr, const char* file, int line, const char* msg = nullptr)
{
    std::stringstream s;
    if (!expr) {
        if (msg) s << msg;
        s << "\nFailed:\t\t" << expr_str << "\n";
        s << "Source:\t\t" << file << ", line " << line << "\n";
        return s.str();
    }
    return {};
}

#define TEST_ASSERT_1(x) __test_assert(#x, x, __FILE__, __LINE__, nullptr)
#define TEST_ASSERT_2(x, msg) __test_assert(#x, x, __FILE__, __LINE__, msg)
#define GET_TEST_ASSERT(_1, _2, NAME, ...) NAME
#define test_assert(...) \
do {\
auto __test_result = GET_TEST_ASSERT(__VA_ARGS__, TEST_ASSERT_2, TEST_ASSERT_1) (__VA_ARGS__); \
if (__test_result) return __test_result; \
} while(0)

#else
#define TEST(test_name, test_body)
#endif
