#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>

#include <fastpm/libfastpm.h>
#include <fastpm/string.h>

char *
fastpm_file_get_content(const char * filename)
{
    FILE * fp = fopen(filename, "r");
    if(!fp) return NULL;

    fseek(fp, 0, SEEK_END);
    size_t file_length = ftell(fp);

    char * buf = malloc(file_length + 1);
    fseek(fp, 0, SEEK_SET);
    fread(buf, 1, file_length, fp);
    fclose(fp);
    buf[file_length] = 0;
    return buf;
}

char **
fastpm_strsplit(const char * str, const char * split)
{
    size_t N = 0;
    const char * p1;
    for(p1 = str; *p1; p1 ++) {
        if(strchr(split, *p1)) N++;
    }
    N++;

    char ** buf = malloc((N + 1) * sizeof(char*) + strlen(str) + 1);
    /* The first part of the buffer is the pointer to the lines */
    /* The second part of the buffer is the actually lines */
    char * dup = (void*) (buf + (N + 1));
    strcpy(dup, str);
    ptrdiff_t i = 0;
    char *p, *q = dup;
    for(p = dup; *p; p ++) {
        if(strchr(split, *p)) {
            buf[i] = q;
            i ++;
            *p = 0;
            q = p + 1;
        }
    }
    buf[i] = q;
    i++;
    buf[i] = NULL;
    return buf;
}

char *
fastpm_strdup(const char * str)
{
    size_t N = strlen(str);
    char * d = malloc(N + 1);
    strcpy(d, str);
    return d;
}

char *
fastpm_strdup_printf(const char * fmt, ...)
{
    va_list va;
    va_start(va, fmt);
    char * buf = fastpm_strdup_vprintf(fmt, va);
    va_end(va);
    return buf;
}

char *
fastpm_strdup_vprintf(const char * fmt, va_list va)
{
    va_list va2;
    va_copy(va2, va);
    /* This relies on a good LIBC vsprintf that returns the number of char */
    char buf0[128];
    size_t N = vsnprintf(buf0, 1, fmt, va);

    char * buf = malloc(N + 100);
    vsnprintf(buf, N + 1, fmt, va2);
    buf[N + 1] = 0;
    va_end(va2);
    return buf;
}
