#ifndef _BIGFILE_H_
#define _BIGFILE_H_
#include <stddef.h>

typedef struct BigFile {
    char * basename;
} BigFile;

typedef struct BigBlockAttr {
    int nmemb;
    char dtype[8];
    char * name;
    char * data;
} BigBlockAttr;

typedef struct BigBlockAttrSet {
    int dirty;
    char * attrbuf;
    size_t bufused;
    size_t bufsize;

    BigBlockAttr * attrlist;
    size_t listused;
    size_t listsize;
} BigBlockAttrSet;

typedef struct BigBlock {
    char dtype[8]; /* numpy style 
                        dtype[0] <: little endian, > big endian, = machine endian
                        dtype[1] type char
                        dtype[2:] width in bytes
                    */
    int nmemb;  /* num of dtype typed elements per item */
    char * basename;
    size_t size;
    size_t * fsize; /* Nfile + 1, in units of elements */
    size_t * foffset; /* Nfile + 1, in units of elements */
    unsigned int * fchecksum; /* sysv sums of each file (unreduced) */
    int Nfile;
    BigBlockAttrSet attrset;
    int dirty;
} BigBlock;

typedef struct BigBlockPtr{
    int fileid;
    ptrdiff_t roffset; /* offset within a file */
    ptrdiff_t aoffset; /* abs offset */ 
} BigBlockPtr;

/**
 *
 * dtype stuff.
 * similar to numpy dtype descr.
 *
 * dtype[0] is endianness:
 *    < LE
 *    > BE
 *    = native
 * dtype[1] is kind
 *    'i'  int
 *    'f'  float
 *    'u'  unsigned int
 * dtype[2:] is the bytewidth
 *
 * dtype[0] can be omitted, in which case native is assumed.
 *
 */

typedef struct BigArray {
    int ndim;
    char dtype[8];
    ptrdiff_t dims[32];
    ptrdiff_t strides[32];
    size_t size;
    void * data;
} BigArray;

typedef struct BigArrayIter {
    ptrdiff_t pos[32];
    BigArray * array;
    int contiguous;
    void * dataptr;
} BigArrayIter;

int big_file_set_buffer_size(size_t bytes);
char * big_file_get_error_message(void);
void big_file_set_error_message(char * msg);
void big_file_checksum(unsigned int * sum, void * buf, size_t size);

int big_file_open(BigFile * bf, const char * basename);
int big_file_create(BigFile * bf, const char * basename);
int big_file_list(BigFile * bf, char *** blocknames, int * Nblocks);
int big_file_open_block(BigFile * bf, BigBlock * block, const char * blockname);
int big_file_create_block(BigFile * bf, BigBlock * block, const char * blockname, const char * dtype, int nmemb, int Nfile, const size_t fsize[]);
int big_file_close(BigFile * bf);
int big_block_flush(BigBlock * block);
int big_file_mksubdir_r(const char * pathname, const char * subdir);

int big_block_open(BigBlock * bb, const char * basename);
int big_block_clear_checksum(BigBlock * bb);
int big_block_create(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, const size_t fsize[]);
int big_block_close(BigBlock * block);
int big_block_seek(BigBlock * bb, BigBlockPtr * ptr, ptrdiff_t offset);
int big_block_seek_rel(BigBlock * bb, BigBlockPtr * ptr, ptrdiff_t rel);
int big_block_eof(BigBlock * bb, BigBlockPtr * ptr);
int big_block_read(BigBlock * bb, BigBlockPtr * ptr, BigArray * array);
int big_block_read_simple(BigBlock * bb, ptrdiff_t start, ptrdiff_t size, BigArray * array, const char * dtype);

int big_block_write(BigBlock * bb, BigBlockPtr * ptr, BigArray * array);
int big_block_remove_attr(BigBlock * block, const char * attrname);
int big_block_set_attr(BigBlock * block, const char * attrname, const void * data, const char * dtype, int nmemb);
int big_block_get_attr(BigBlock * block, const char * attrname, void * data, const char * dtype, int nmemb);
BigBlockAttr * big_block_lookup_attr(BigBlock * block, const char * attrname);
BigBlockAttr * big_block_list_attrs(BigBlock * block, size_t * count);

int dtype_normalize(char * dst, const char * src);
void dtype_format(char * buffer, const char * dtype, const void * data, const char * flags);
void dtype_parse(const char * buffer, const char * dtype, void * data, const char * fmt);
int dtype_convert(BigArrayIter * dst, BigArrayIter * src, size_t nmemb);
int dtype_convert_simple(void * dst, const char * dstdtype, const void * src, const char * srcdtype, size_t nmemb);
int dtype_cmp(const char * dtype1, const char * dtype2);
char dtype_kind(const char * dtype);
int dtype_needswap(const char * dtype);
int dtype_itemsize(const char * dtype);

int big_array_init(BigArray * array, void * buf, const char * dtype, int ndim, const size_t dims[], const ptrdiff_t strides[]);
int big_array_iter_init(BigArrayIter * iter, BigArray * array);
void big_array_iter_advance(BigArrayIter * iter);
#endif
