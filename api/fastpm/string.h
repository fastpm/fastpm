#include <stdarg.h>

char *
fastpm_file_get_content(const char * filename);

char **
fastpm_strsplit(const char * str, const char * split);

char *
fastpm_strdup(const char * str);

char *
fastpm_strdup_printf(const char * fmt, ...);

char *
fastpm_strdup_vprintf(const char * fmt, va_list va);
