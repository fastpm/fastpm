#ifndef _FASTPM_STRING_H_
#define _FASTPM_STRING_H_

#include <stdarg.h>

FASTPM_BEGIN_DECLS

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

void
fastpm_path_ensure_dirname(const char * path);

FASTPM_END_DECLS
#endif
