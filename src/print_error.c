#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include "print_error.h"

void print_error(bool use_errno, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);

    fprintf(stderr, "[ERROR] ");

    vfprintf(stderr, fmt, args);

    if (use_errno) {
        char buf[256];
#if defined(_GNU_SOURCE) && !defined(__APPLE__)
        char *err_str = strerror_r(errno, buf, sizeof(buf));
        fprintf(stderr, ": %s", err_str);
#else
        if (strerror_r(errno, buf, sizeof(buf)) == 0) {
            fprintf(stderr, ": %s", buf);
        } else {
            fprintf(stderr, ": Unknown error %d", errno);
        }
#endif
    }

    fprintf(stderr, "\n");
    va_end(args);
}
