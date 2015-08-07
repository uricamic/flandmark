#ifndef FLANDMAR_MSVC_COMPAT
#define FLANDMAR_MSVC_COMPAT
#ifdef _MSC_VER

typedef unsigned char uint8_t;
typedef char int8_t;
typedef unsigned __int16 uint16_t;
typedef __int16 int16_t;
typedef unsigned __int32 uint32_t;
typedef __int32 int32_t;
typedef unsigned __int64 uint64_t;
typedef __int64 int64_t;

#else
#	include <stdint.h>
#endif

#endif
