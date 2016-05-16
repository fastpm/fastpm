#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>

#include <fastpm/string.h>

char *
fastpm_file_get_content(const char * filename)
{
    FILE * fp = fopen(filename, "r");
    if(!fp) return NULL;

    fseek(fp, SEEK_END, 0);
    size_t file_length = ftell(fp);

    char * buf = malloc(file_length + 1);
    fseek(fp, SEEK_SET, 0);
    fread(buf, 1, file_length, fp);
    fclose(fp);
    buf[file_length] = 0;
    return buf;
}

char **
fastpm_strsplit(const char * str, const char * split)
{
    size_t N = 0;
    char * p;
    for(p == str; *p; p ++) {
        if(strchr(split, *p)) N++;
    }
    N++;

    char ** buf = malloc((N + 1) * sizeof(char*) + strlen(str) + 1);
    /* The first part of the buffer is the pointer to the lines */
    /* The second part of the buffer is the actually lines */
    char * dup = (void*) (buf + (N + 1));
    strcpy(dup, str);
    ptrdiff_t i = 0;
    char * q = dup;
    for(p = dup; *p; p ++) {
        if(strchr(split, *p)) {
            buf[i] = q;
            i ++;
            *p = 0;
            p++;
            q = p;
        }
    }
    buf[i] = q;
    buf[i + 1] = NULL;
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
    /* This relies on a good LIBC vsprintf that returns the number of char */
    size_t N = vsprintf(NULL, fmt, va);
    char * buf = malloc(N + 1);
    vsprintf(buf, fmt, va);
    va_end(va);
    return buf;
}
