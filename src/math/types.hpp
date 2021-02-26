#include <cmath>

#ifdef __GNUC__
#include <quadmath.h>
#endif

namespace dptk {

typedef float  b32;
typedef double b64;

#ifdef __GNUC__
typedef __float128 b128;
#endif

typedef bool               u1;
typedef unsigned char      u8;
typedef unsigned short int u16;
typedef unsigned int       u32;
typedef unsigned long int  u64;

#ifdef __GNUC__
typedef __uint128_t u128;
#endif

typedef char      i8;
typedef short int i16;
typedef int       i32;
typedef long int  i64;

#ifdef __GNUC__
typedef __int128_t i128;
#endif

} // namespace dptk
