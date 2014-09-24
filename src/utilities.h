#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <execinfo.h>

#include <stdarg.h>
#include <string.h>
#include <time.h>

#include <errno.h>
#include <sys/types.h> 
#include <zlib.h>

////////////////////////////////////////////////////////////////////////////////
// These are just wrappers around the C assert. forceassert is a way to enforce
// asserts in the production build as well. Only inexpensive and vital asserts
// should be encoded with forceassert
////////////////////////////////////////////////////////////////////////////////

#define pre  assert
#define post assert

#define ForceAssertWithDebug(e, f, l)    \
    if (!(e))                        \
    {                                \
        fprintf(stderr, "Assertion failed: %s file %s line %d\n", #e, f, l);\
        exit(EXIT_FAILURE);          \
    }

#define ForceAssert(expr) ForceAssertWithDebug(expr, __FILE__, __LINE__)

// end of asserts block

////////////////////////////////////////////////////////////////////////////////
// constants defined throughout
////////////////////////////////////////////////////////////////////////////////

// bool constants
#ifndef __cplusplus
    typedef enum Bool_em {
        FALSE = 0,
        TRUE  = 1
    }Bool;
#else
    typedef bool Bool;
    #define TRUE true
    #define FALSE false
#endif

// some systems do not have uint, ulong, uchar and uint128_t
typedef unsigned int uint;
typedef unsigned char uchar;
typedef __uint128_t uint128_t;

// valid nucleotides
#define _ -1

static const signed char fasta_encoding[] = {
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,0,_,1,_,_,
_,2,_,_,_,_,_,_,0,_,
_,_,_,_,3,_,_,_,_,_,
_,_,_,_,_,_,_,0,_,1,
_,_,_,2,_,_,_,_,_,_,
0,_,_,_,_,_,3,_,_,_,
_,_,_,_,_,_,_,_,_
};

static const char bit_encoding[] = {'A', 'C', 'G', 'T'};

// ignore warnings for these local variables
# define UNUSED __attribute__((unused))

// maximum and minimum values out of two
#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// maximum value of a particular integer datatype
#define umaxof(t) (((0x1ULL << ((sizeof(t) * 8ULL) - 1ULL)) - 1ULL) | \
                    (0xFULL << ((sizeof(t) * 8ULL) - 4ULL)))

#define smaxof(t) (((0x1ULL << ((sizeof(t) * 8ULL) - 1ULL)) - 1ULL) | \
                    (0x7ULL << ((sizeof(t) * 8ULL) - 4ULL)))

// end of constants

////////////////////////////////////////////////////////////////////////////////
// Error reporting routines
////////////////////////////////////////////////////////////////////////////////

extern time_t t0;  // time when we started the execution of this program

void _fatalf(const char* const filename,
             const int lineNumber,
             const char* const fmt, ...);
void _warnf(const char* const filename,
             const int lineNumber,
             const char* const fmt, ...);
void _debugf(const char* const filename,
             const int lineNumber,
             const char* const fmt, ...);
void _fatal(const char* const filename,
            const int lineNumber,
            const char* const msg);

#define PrintThenDie(x) _fatal(__FILE__, __LINE__, x)
#define PrintMessageThenDie(x, ...) _fatalf(__FILE__, __LINE__, x, __VA_ARGS__)
#define PrintWarning(x, ...) _warnf(__FILE__, __LINE__, x, __VA_ARGS__)
#define PrintDebugMessage(x, ...) _debugf(__FILE__, __LINE__, x, __VA_ARGS__)

// end of error reporting routines

////////////////////////////////////////////////////////////////////////////////
// memory allocation routines
////////////////////////////////////////////////////////////////////////////////

void* _ckalloc(const size_t size,
               const char* const filename,
               const uint linenum);
void* _ckallocz(const size_t size,
                const char* const filename,
                const uint linenum);
void* _ckrealloc(void* ptr,
                 const size_t size,
                 const char* const filename,
                 const uint linenum);
void _ckfree(void* ptr,
             const char* const filename,
             const uint linenum);

// Routines for memory allocation
#define CkallocOrDie(x) _ckalloc(x, __FILE__, __LINE__)
#define CkalloczOrDie(x) _ckallocz(x, __FILE__, __LINE__)
#define CkreallocOrDie(x, y) _ckrealloc(x, y, __FILE__, __LINE__)

// Routine to deallocate memory
#define Ckfree(x) _ckfree(x, __FILE__, __LINE__)

// end of memory allocation routines

////////////////////////////////////////////////////////////////////////////////
// file handling routines
////////////////////////////////////////////////////////////////////////////////

FILE* CkopenOrDie(const char* const name, const char* const mode);

// a drop in replacement for getline
ssize_t Getline(char** lineptr, size_t* max, FILE* stream);

// a function similar to getline, but for zipped files
ssize_t GetZippedLine(char** lineptr, size_t* max, gzFile* stream);

// end of file handling routines

////////////////////////////////////////////////////////////////////////////////
// string handling routines
////////////////////////////////////////////////////////////////////////////////

// return a copy of this string 
char* CopyString(const char* const str);

// determine whether two strings are identical 
Bool SameString(const char *s, const char *t);

// determine whether the first len letters of the two strings are the same? 
Bool CompareNames(const char* const s1, const char* const s2, const int len);

// reverse complement the string in place 
char* ReverseComplementString(char* sequence, const int length);

// reverse string in place 
char* ReverseString(char* sequence, const int length);

// end of string handling routines

////////////////////////////////////////////////////////////////////////////////
// miscellaneous routines          
////////////////////////////////////////////////////////////////////////////////

// this is a routine to report the rough memory usage of the program at the
// current time point. 
void _report(const char* const filename,
             const int lineNumber);

#define ReportMemoryUsage() _report(__FILE__, __LINE__)

#endif
