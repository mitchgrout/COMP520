#ifndef __UTILS
#define __UTILS
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

// Log data to stream
static inline void log(FILE *stream, const char *fmt, ...)
{
    // Print timestamp
    {
        char buf[8+1];
        time_t t;
        time(&t);
        struct tm *info = localtime(&t);
        strftime(buf, 9, "%H:%M:%S", info);
        fprintf(stream, "[%s] ", buf);
    }
    // Print fmt
    {
        va_list args;
        va_start(args, fmt);
        vfprintf(stream, fmt, args);
        fputc('\n', stream);
        va_end(args);
    }
}
#endif // __UTILS
