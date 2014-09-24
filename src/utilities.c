#include "utilities.h"

time_t t0;

void _fatalf(const char* const filename,
             const int lineNumber,
             const char* const fmt, ...) {
    pre(fmt != NULL);

    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    fprintf(stderr, "%s:%d ", filename, lineNumber);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(EXIT_FAILURE);
}

void _warnf(const char* const filename,
            const int lineNumber,
            const char* const fmt, ...) {
    pre(fmt != NULL);

    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    fprintf(stderr, "%s:%d ", filename, lineNumber);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
}

void _debugf(const char* const filename __attribute__((unused)),
             const int lineNumber __attribute__((unused)),
             const char* const fmt, ...) {
    pre(fmt != NULL);

    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    time_t t1 = time(0);
    fprintf(stderr, "[%.0f s]\t", difftime(t1,t0));
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
}

/*print the message on stderr and die*/
void _fatal(const char* const filename,
            const int lineNumber,
            const char* const msg) {
    pre(msg != NULL);

    _fatalf(filename, lineNumber, "%s", msg);
}

// make a call to malloc, and panic if the allocation request fails.
void* _ckalloc(const size_t size,
              const char* const filename UNUSED,
              const uint linenum UNUSED) {
    void* ptr;
    if ((ptr = malloc(size)) == NULL) {
        fprintf(stderr, "Error in allocating %zd bytes.\n", size);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

#ifdef DebugMemory
    fprintf(stderr,
    "%s:%d Allocated %zd bytes at %p.\n", filename, linenum, size, ptr);
#endif

    return ptr;
}

// fill the allocated memory with '0' before returning a pointer to it
void* _ckallocz(const size_t size,
               const char* const filename UNUSED,
               const uint linenum UNUSED) {
    void* ptr = _ckalloc(size, filename, linenum);

    if (memset(ptr, 0, size) != ptr) {
        fprintf(stderr, "Error in initializing the area\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

/*reallocate the memory to the given size*/
void* _ckrealloc(void* p,
                 const size_t size,
                 const char* const filename UNUSED,
                 const uint linenum UNUSED) {
    void* ptr UNUSED = p;
    p = p ? realloc(p, size) : malloc(size);
    if (!p) {
        PrintThenDie("ckrealloc failed");
    }

#ifdef DebugMemory
    fprintf(stderr,
    "%s:%d Reallocated %zd bytes from %p to %p\n",
    filename, linenum, size, ptr, p);
#endif

    return p;
}

// free the memory pointed to by this pointer.
void _ckfree(void* ptr,
             const char* const filename UNUSED,
             const uint linenum UNUSED) {
#ifdef DebugMemory
    fprintf(stderr, "%s:%d Deallocated bytes at %p.\n", filename, linenum, ptr);
#endif

    if (ptr) free(ptr);
}


FILE* CkopenOrDie(const char* const name, const char* const mode) {
    FILE* fp;

    if ((fp = fopen(name, mode)) == NULL) {
        PrintMessageThenDie("error in opening the file %s: %s",
        name, strerror(errno));
    }
    return fp;
}

ssize_t Getline(char** lineptr, size_t* max, FILE* stream) {
    int ch;
    ssize_t size = 0;
    char* ptr;
    size_t sz;
    
    while ((ch = fgetc(stream)) != EOF) {
        if (size >= *(ssize_t*)max) {
            sz = size + (size >> 5) + 16;
            *max = sz;
            if ((ptr = CkreallocOrDie(*lineptr, *max)) == NULL) {
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size++] = ch;
        if (ch == '\n') {
            break;
        }
    }
    
    if (size != 0) {
        if(size >= *(ssize_t*)max){
            sz = size + (size >> 5) + 16;
            *max = sz;
            if ((ptr = CkreallocOrDie(*lineptr, *max)) == NULL) {
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size] = '\0';
    }
    if (0 == size || ch == EOF) {
        return -1;
    }
    
    return size;
}

// a function similar to getline, but for zipped files
ssize_t GetZippedLine(char** lineptr, size_t* max, gzFile* stream) {
    int ch;
    ssize_t size = 0;
    char* ptr;
    size_t sz;

    while ((ch = gzgetc(*stream)) != EOF) {
        if (size >= *(ssize_t*)max) {
            sz = size + (size >> 5) + 16;
            *max = sz;
            if ((ptr = CkreallocOrDie(*lineptr, *max)) == NULL) {
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size++] = ch;
        if (ch == '\n') {
            break;
        }
    }
    
    if (size != 0) {
        if(size >= *(ssize_t*)max){
            sz = size + (size >> 5) + 16;
            *max = sz;
            if ((ptr = CkreallocOrDie(*lineptr, *max)) == NULL) {
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size] = '\0';
    }
    if (0 == size || ch == EOF) {
        return -1;
    }

    return size;
}

const unsigned char dna_complement[] =
  "                                                                "
  " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
  "                                                                "
  "                                                                ";

/**-------------------String handling routines------------------------------*/

/* return a copy of the string */
char* CopyString(const char* const str) {
    char* copy = CkalloczOrDie(strlen(str)+1);
    snprintf(copy, strlen(str) + 1, "%s", str);
    return copy;
}

/* determine whether two strings are identical */
Bool SameString(const char *s, const char *t) {
    return (strcmp(s, t) == 0);
}

/* determine whether two strings are identical for the first len characters */
Bool CompareNames(const char* const s1, const char* const s2, const int len) {
    int i;

    for (i = 0; i < len; i++) {
        if (*(s1+i) != *(s2+i)) {
            return FALSE;
        }
    }

    return TRUE;
}

/*reverse complement the string in place*/
char* ReverseComplementString(char* sequence, const int length) {
    char* s = sequence;
    char* p = s + length - 1;

    while (s <= p) {
        char c;

        c  = dna_complement[(int)*s];
        *s = dna_complement[(int)*p];
        *p = c;
        ++s; --p;
    }
    return sequence;
}

/*reverse string in place*/
char* ReverseString(char* sequence, const int length) {
    char* s = sequence;
    char* p = s + length - 1;

    while (s <= p) {
        char c;

         c = *s;
        *s = *p;
        *p =  c;
        ++s; --p;
    }
    return sequence;
}

/* report the memory usage of the program at this instant */
void _report(const char* const filename,
             const int lineNumber) {
    char cmd[1024];
    sprintf(cmd, "ps up %d | awk '{print $6}'", getpid());
    FILE* fp = popen(cmd, "r");
    if (fp == NULL) { 
        // popen fails, then I just return
        return;
    }
    char path[1024];
    while (fgets(path, sizeof(path)-1, fp) != NULL) {
    }
    pclose(fp);
    path[strlen(path)] = 0;

    fprintf(stderr, "[%s:%d] Memory used (in kB) : %s", 
    filename, lineNumber, path);
}
