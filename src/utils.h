#pragma once

#include <iostream>
#include <string>

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

typedef f32 scalar;
//TO-DO define what matrix type we use globaly
//typedef eingen::matrix matrix

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

#define GB_DEFER_1(x, y) x##y
#define GB_DEFER_2(x, y) GB_DEFER_1(x, y)
#define GB_DEFER_3(x)    GB_DEFER_2(x, __COUNTER__)
#define defer(code)      auto GB_DEFER_3(_defer_) = gb__defer_func([&]()->void{code;})

#ifdef NDEBUG
#define assert(...)
#else
#define GB_ASSERT_1(x) gb__assert(#x, x, __FILE__, __LINE__, nullptr)
#define GB_ASSERT_2(x, msg) gb__assert(#x, x, __FILE__, __LINE__, msg)
#define GET_GB_ASSERT(_1, _2, NAME, ...) NAME
#define assert(...) GET_GB_ASSERT(__VA_ARGS__, GB_ASSERT_2, GB_ASSERT_1) (__VA_ARGS__)

inline void gb__assert(const char* expr_str, bool expr, const char* file, int line, const char* msg = nullptr)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t";
        if (msg) std::cerr << msg;
        std::cerr << "\nExpected:\t" << expr_str << "\n" << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}
#endif
