#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

#include "bigfile.h"

#define EXT_HEADER "header"
#define EXT_ATTR "attr"
#define EXT_ATTR_V2 "attr-v2"
#define EXT_DATA   "%06X"
#define FILEID_ATTR -2
#define FILEID_ATTR_V2 -3
#define FILEID_HEADER -1

#define RAISE(ex, errormsg, ...) { __raise__(errormsg, __FILE__, __LINE__, ##__VA_ARGS__); goto ex; } 
#define RAISEIF(condition, ex, errormsg, ...) { if(condition) RAISE(ex, errormsg, ##__VA_ARGS__); }

static char * ERRORSTR = NULL;
static size_t CHUNK_BYTES = 64 * 1024 * 1024;

/* Internal AttrSet API */

struct BigAttrSet {
    /* All members are readonly */
    int dirty;
    char * attrbuf;
    size_t bufused;
    size_t bufsize;

    BigAttr * attrlist;
    size_t listused;
    size_t listsize;
};

static BigAttrSet *
attrset_create();
static void attrset_free(BigAttrSet * attrset);
static int
attrset_read_attr_set_v1(BigAttrSet * attrset, const char * basename);
static int
attrset_read_attr_set_v2(BigAttrSet * attrset, const char * basename);
static int
attrset_write_attr_set_v2(BigAttrSet * attrset, const char * basename);
static BigAttr *
attrset_lookup_attr(BigAttrSet * attrset, const char * attrname);
static int
attrset_remove_attr(BigAttrSet * attrset, const char * attrname);
static BigAttr *
attrset_list_attrs(BigAttrSet * attrset, size_t * count);
static int
attrset_set_attr(BigAttrSet * attrset, const char * attrname, const void * data, const char * dtype, int nmemb);
static int
attrset_get_attr(BigAttrSet * attrset, const char * attrname, void * data, const char * dtype, int nmemb);

/* Internal dtype API */
static int
dtype_convert_simple(void * dst, const char * dstdtype, const void * src, const char * srcdtype, size_t nmemb);

/* global settings */
int
big_file_set_buffer_size(size_t bytes)
{
    CHUNK_BYTES = bytes;
    return 0;
}

/* Error handling */
char * big_file_get_error_message() {
    return ERRORSTR;
}
void big_file_clear_error_message() {
    if(ERRORSTR) {
        free(ERRORSTR);
        ERRORSTR = NULL;
    }
}

static char * _path_join(const char * part1, const char * part2)
{
    size_t l1 = strlen(part1);
    size_t l2 = strlen(part2);
    char * result = malloc(l1 + l2 + 10);
    if(l1 == 0) {
        sprintf(result, "%s", part2);
    } else
    if(part1[l1 - 1] == '/') {
        sprintf(result, "%s%s", part1, part2);
    } else {
        sprintf(result, "%s/%s", part1, part2);
    }
    return result;
}

static int
_big_file_path_is_block(const char * basename)
{
    FILE * fheader = _big_file_open_a_file(basename, FILEID_HEADER, "r");
    if(fheader == NULL) return 0;
    fclose(fheader);
    return 1;
}


static void 
__raise__(const char * msg, const char * file, const int line, ...)
{
    char * mymsg;
    if(!msg) {
        if(ERRORSTR) {
            mymsg = strdup(ERRORSTR);
        } else {
            mymsg = strdup("UNKNOWN ERROR");
        }
    } else {
        mymsg = strdup(msg);
    }

    if(ERRORSTR) free(ERRORSTR);
    ERRORSTR = malloc(strlen(mymsg) + 512);

    char * fmtstr = alloca(strlen(mymsg) + 512);
    sprintf(fmtstr, "%s @(%s:%d)", mymsg, file, line);
    free(mymsg);
    va_list va;
    va_start(va, line);
    vsprintf(ERRORSTR, fmtstr, va); 
    va_end(va);
}

void
big_file_set_error_message(char * msg)
{
    if(ERRORSTR) free(ERRORSTR);
    ERRORSTR = strdup(msg);
}

/* BigFile */

int
big_file_open(BigFile * bf, const char * basename)
{
    memset(bf, 0, sizeof(bf[0]));
    struct stat st;
    RAISEIF(0 != stat(basename, &st),
            ex_stat,
            "Big File `%s' does not exist (%s)", basename,
            strerror(errno));
    bf->basename = strdup(basename);
    return 0;
ex_stat:
    return -1;
}

int big_file_create(BigFile * bf, const char * basename) {
    memset(bf, 0, sizeof(bf[0]));
    bf->basename = strdup(basename);
    RAISEIF(0 != _big_file_mksubdir_r("", basename),
        ex_subdir,
        NULL);
    return 0;
ex_subdir:
    return -1;
}

struct bblist {
    struct bblist * next;
    char * blockname;
};
static int (filter)(const struct dirent * ent) {
    if(ent->d_name[0] == '.') return 0;
    return 1;
}
static struct
bblist * listbigfile_r(const char * basename, char * blockname, struct bblist * bblist) {
    struct dirent **namelist;
    struct stat st;
    int n;
    int i;

    char * current;
    current = _path_join(basename, blockname);

    n = scandir(current, &namelist, filter, alphasort);
    free(current);
    if (n < 0)
        return bblist;

    for(i = 0; i < n ; i ++) {
        char * blockname1 = _path_join(blockname, namelist[i]->d_name);
        char * fullpath1 = _path_join(basename, blockname1);
        BigBlock bb = {0};
        if(_big_file_path_is_block(fullpath1)) {
            struct bblist * n = malloc(sizeof(struct bblist) + strlen(blockname1) + 1);
            n->next = bblist;
            n->blockname = (char*) &n[1];
            strcpy(n->blockname, blockname1);
            bblist = n;
        } else {
            big_file_clear_error_message();
            bblist = listbigfile_r(basename, blockname1, bblist);
        }
        free(fullpath1);
        free(blockname1);
        free(namelist[i]);
    }
    free(namelist);
    return bblist;
}

int
big_file_list(BigFile * bf, char *** blocknames, int * Nblocks)
{
    struct bblist * bblist = listbigfile_r(bf->basename, "", NULL);
    struct bblist * p;
    int N = 0;
    int i;
    for(p = bblist; p; p = p->next) {
        N ++;
    }
    *Nblocks = N;
    *blocknames = malloc(sizeof(char*) * N);
    for(i = 0; i < N; i ++) {
        (*blocknames)[i] = strdup(bblist->blockname);
        p = bblist;
        bblist = bblist->next;
        free(p); 
    }
    return 0;
}

int
big_file_open_block(BigFile * bf, BigBlock * block, const char * blockname)
{
    char * basename = _path_join(bf->basename, blockname);
    int rt = big_block_open(block, basename);
    free(basename);
    return rt;
}

int
big_file_create_block(BigFile * bf, BigBlock * block, const char * blockname, const char * dtype, int nmemb, int Nfile, const size_t fsize[])
{
    RAISEIF(0 != _big_file_mksubdir_r(bf->basename, blockname),
            ex_subdir,
            NULL);
    char * basename = _path_join(bf->basename, blockname);
    int rt = big_block_create(block, basename, dtype, nmemb, Nfile, fsize);
    free(basename);
    return rt;
ex_subdir:
    return -1;
}

int
big_file_close(BigFile * bf)
{
    free(bf->basename);
    bf->basename = NULL;
    return 0;
}
static void
sysvsum(unsigned int * sum, void * buf, size_t size);

void
big_file_checksum(unsigned int * sum, void * buf, size_t size)
{
    sysvsum(sum, buf, size);
}


/* Bigblock */

int
big_block_open(BigBlock * bb, const char * basename)
{
    memset(bb, 0, sizeof(bb[0]));
    if(basename == NULL) basename = "";
    bb->basename = strdup(basename);
    bb->dirty = 0;
    FILE * fheader = _big_file_open_a_file(bb->basename, FILEID_HEADER, "r");
    RAISEIF (!fheader,
            ex_open,
            NULL);

    RAISEIF(
           (1 != fscanf(fheader, " DTYPE: %s", bb->dtype)) ||
           (1 != fscanf(fheader, " NMEMB: %d", &(bb->nmemb))) ||
           (1 != fscanf(fheader, " NFILE: %d", &(bb->Nfile))),
           ex_fscanf,
           "Failed to read header of block `%s' (%s)", bb->basename, strerror(errno));

    bb->fsize = calloc(bb->Nfile, sizeof(size_t));
    RAISEIF(!bb->fsize,
            ex_fsize,
            "Failed to alloc memory of block `%s'", bb->basename);
    bb->foffset = calloc(bb->Nfile + 1, sizeof(size_t));
    RAISEIF(!bb->foffset,
            ex_foffset,
            "Failed to alloc memory of block `%s'", bb->basename);
    bb->fchecksum = calloc(bb->Nfile, sizeof(int));
    RAISEIF(!bb->fchecksum,
            ex_fchecksum,
            "Failed to alloc memory `%s'", bb->basename);
    int i;
    for(i = 0; i < bb->Nfile; i ++) {
        int fid; 
        size_t size;
        unsigned int cksum;
        unsigned int sysv;
        RAISEIF(4 != fscanf(fheader, " " EXT_DATA ": %td : %u : %u", &fid, &size, &cksum, &sysv),
                ex_fscanf1,
                "Failed to readin physical file layout `%s' %d (%s)", bb->basename, fid,
                strerror(errno));
        bb->fsize[fid] = size;
        bb->fchecksum[fid] = cksum;
    }
    bb->foffset[0] = 0;
    for(i = 0; i < bb->Nfile; i ++) {
        bb->foffset[i + 1] = bb->foffset[i] + bb->fsize[i];
    }
    bb->size = bb->foffset[bb->Nfile];

    bb->attrset = attrset_create();

    /* COMPATIBLE WITH V1 ATTR FILES */
    RAISEIF(0 != attrset_read_attr_set_v1(bb->attrset, bb->basename),
            ex_readattr,
            NULL);

    RAISEIF(0 != attrset_read_attr_set_v2(bb->attrset, bb->basename),
            ex_readattr,
            NULL);

    fclose(fheader);
    return 0;

ex_readattr:
ex_fscanf1:
    free(bb->fchecksum);
ex_fchecksum:
    free(bb->foffset);
ex_foffset:
    free(bb->fsize);
ex_fsize:
ex_fscanf:
    fclose(fheader);
ex_open:
    return -1;
}

int
_big_block_create_internal(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, const size_t fsize[])
{
    memset(bb, 0, sizeof(bb[0]));
    if(basename == NULL) basename = "";
    bb->basename = strdup(basename);

    if(dtype == NULL) {
        dtype = "i8";
        Nfile = 0;
        fsize = NULL;
    }
    /* always use normalized dtype in files. */
    dtype_normalize(bb->dtype, dtype);

    bb->Nfile = Nfile;
    bb->nmemb = nmemb;
    bb->fsize = calloc(bb->Nfile, sizeof(size_t));
    RAISEIF(!bb->fsize, ex_fsize, "No memory"); 
    bb->fchecksum = calloc(bb->Nfile, sizeof(int));
    RAISEIF(!bb->fchecksum, ex_fchecksum, "No memory"); 
    bb->foffset = calloc(bb->Nfile + 1, sizeof(size_t));
    RAISEIF(!bb->foffset, ex_foffset, "No memory"); 
    int i;
    bb->foffset[0] = 0;
    for(i = 0; i < bb->Nfile; i ++) {
        bb->fsize[i] = fsize[i];
        bb->foffset[i + 1] = bb->foffset[i] + bb->fsize[i];
        bb->fchecksum[i] = 0;
    }

    bb->size = bb->foffset[bb->Nfile];

    bb->attrset = attrset_create();
    bb->attrset->dirty = 1;
    bb->dirty = 1;
    RAISEIF(0 != big_block_flush(bb), 
            ex_flush, NULL);

    bb->dirty = 0;

    return 0;

ex_flush:
ex_fileio:
    free(bb->foffset);
ex_foffset:
    free(bb->fchecksum);
ex_fchecksum:
    free(bb->fsize);
ex_fsize:
    return -1;
}

static void
big_block_close_internal(BigBlock * block);

int
big_block_create(BigBlock * bb, const char * basename, const char * dtype, int nmemb, int Nfile, const size_t fsize[])
{
    int rt = _big_block_create_internal(bb, basename, dtype, nmemb, Nfile, fsize);
    int i;
    RAISEIF(rt != 0,
                ex_internal,
                NULL);

    /* now truncate all files */
    for(i = 0; i < bb->Nfile; i ++) {
        FILE * fp = _big_file_open_a_file(bb->basename, i, "w");
        RAISEIF(fp == NULL, 
                ex_fileio, 
                NULL);
        fclose(fp);
    }
ex_internal:
    return rt;
ex_fileio:
    big_block_close_internal(bb);
    return -1;
}

int
big_block_clear_checksum(BigBlock * bb)
{
    memset(bb->fchecksum, 0, bb->Nfile * sizeof(int));
    return 0;
}

void
big_block_set_dirty(BigBlock * block, int value)
{
    block->dirty = value;
}

int
big_block_flush(BigBlock * block)
{
    FILE * fheader = NULL;
    if(block->dirty) {
        int i;
        fheader = _big_file_open_a_file(block->basename, FILEID_HEADER, "w+");
        RAISEIF(fheader == NULL, ex_fileio, NULL);
        RAISEIF(
            (0 > fprintf(fheader, "DTYPE: %s\n", block->dtype)) ||
            (0 > fprintf(fheader, "NMEMB: %d\n", block->nmemb)) ||
            (0 > fprintf(fheader, "NFILE: %d\n", block->Nfile)),
                ex_fprintf,
                "Writing file header");
        for(i = 0; i < block->Nfile; i ++) {
            unsigned int s = block->fchecksum[i];
            unsigned int r = (s & 0xffff) + ((s & 0xffffffff) >> 16);
            unsigned int checksum = (r & 0xffff) + (r >> 16);
            RAISEIF(0 > fprintf(fheader, EXT_DATA ": %td : %u : %u\n", i, block->fsize[i], block->fchecksum[i], checksum),
                ex_fprintf, "Writing file information to header");
        }
        fclose(fheader);
        block->dirty = 0;
    }
    if(block->attrset->dirty) {
        RAISEIF(0 != attrset_write_attr_set_v2(block->attrset, block->basename),
            ex_write_attr,
            NULL);
        block->attrset->dirty = 0;
    }
    return 0;

ex_fprintf:
    fclose(fheader);
ex_write_attr:
ex_fileio:
    return -1;
}

static void
big_block_close_internal(BigBlock * block)
{
    attrset_free(block->attrset);

    free(block->basename);
    free(block->fchecksum);
    free(block->fsize);
    free(block->foffset);
    memset(block, 0, sizeof(BigBlock));
}

int
big_block_close(BigBlock * block)
{
    int rt = 0;
    RAISEIF(0 != big_block_flush(block),
            ex_flush,
            NULL);
finalize:
    big_block_close_internal(block);
    return rt;

ex_flush:
    rt = -1;
    goto finalize;
}

BigAttr *
big_block_lookup_attr(BigBlock * block, const char * attrname)
{
    return attrset_lookup_attr(block->attrset, attrname);
}

int
big_block_remove_attr(BigBlock * block, const char * attrname)
{
    return attrset_remove_attr(block->attrset, attrname);
}

BigAttr *
big_block_list_attrs(BigBlock * block, size_t * count)
{
    return attrset_list_attrs(block->attrset, count);
}
int
big_block_set_attr(BigBlock * block, const char * attrname, const void * data, const char * dtype, int nmemb)
{
    return attrset_set_attr(block->attrset, attrname, data, dtype, nmemb);
}

int
big_block_get_attr(BigBlock * block, const char * attrname, void * data, const char * dtype, int nmemb)
{
    return attrset_get_attr(block->attrset, attrname, data, dtype, nmemb);
}

/* *
 * seek ptr to offset, on bb
 *  
 *      offset: 0 : Start
 *              -1 : End
 *              -2 : End - 1
 *
 * returns 0
 *
 * 0 4 5 10 140  
 * */
int
big_block_seek(BigBlock * bb, BigBlockPtr * ptr, ptrdiff_t offset)
{
    /* handle 0 sized files */
    if(bb->size == 0 && offset == 0) {
        ptr->fileid = 0;
        ptr->roffset = 0;
        ptr->aoffset = 0;
        return 0;
    }
    /* handle negatives */
    if(offset < 0) offset += bb->foffset[bb->Nfile];

    RAISEIF(offset > bb->size, 
            ex_eof,
        /* over the end of file */
        /* note that we allow seeking at the end of file */
            "Over the end of file");
    {
        int left = 0;
        int right = bb->Nfile;
        while(right > left + 1) {
            int mid = ((right - left) >> 1) + left;
            if(bb->foffset[mid] <= offset) {
                left = mid;
            } else {
                right = mid;
            }
        }
        ptr->fileid = left;
        ptr->roffset = offset - bb->foffset[left];
        ptr->aoffset = offset;
        return 0;
    }
ex_eof:
    return -1;
}

int
big_block_seek_rel(BigBlock * bb, BigBlockPtr * ptr, ptrdiff_t rel)
{
    ptrdiff_t abs = bb->foffset[ptr->fileid] + ptr->roffset + rel;
    return big_block_seek(bb, ptr, abs);
}

int
big_block_eof(BigBlock * bb, BigBlockPtr * ptr)
{
    ptrdiff_t abs = bb->foffset[ptr->fileid] + ptr->roffset;
    return abs >= bb->size;
}

/* 
 * this function will alloc memory in array and read from offset start 
 * size of rows from the block. 
 * free(array->data) after using it.
 *
 * at most size rows are read, array.dims[0] has the number that's read.
 *
 * if dtype is NULL use the dtype of the block.
 * otherwise cast the array to the dtype
 * */
int
big_block_read_simple(BigBlock * bb, ptrdiff_t start, ptrdiff_t size, BigArray * array, const char * dtype)
{
    BigBlockPtr ptr = {0};
    if(dtype == NULL) {
        dtype = bb->dtype;
    }
    void * buffer;
    size_t dims[2];

    RAISEIF(0 != big_block_seek(bb, &ptr, start),
       ex_seek,
       "failed to seek"       
    );

    if(start + size > bb->size){
        size = bb->size - start;
    }
    RAISEIF(size < 0,
            ex_seek,
            "failed to seek");

    buffer = malloc(size * dtype_itemsize(dtype) * bb->nmemb);

    dims[0] = size;
    dims[1] = bb->nmemb;

    big_array_init(array, buffer, dtype, 2, dims, NULL);

    RAISEIF(0 != big_block_read(bb, &ptr, array),
            ex_read,
            "failed to read");
    return 0;
ex_read:
    free(buffer);
    array->data = NULL;
ex_seek:
    return -1;
}

int
big_block_read(BigBlock * bb, BigBlockPtr * ptr, BigArray * array)
{
    char * chunkbuf = malloc(CHUNK_BYTES);
    int felsize = dtype_itemsize(bb->dtype) * bb->nmemb;
    size_t CHUNK_SIZE = CHUNK_BYTES / felsize;

    BigArray chunk_array = {0};
    size_t dims[2];
    dims[0] = CHUNK_SIZE;
    dims[1] = bb->nmemb;

    BigArrayIter chunk_iter;
    BigArrayIter array_iter;

    FILE * fp = NULL;
    ptrdiff_t toread = 0;

    RAISEIF(chunkbuf == NULL,
            ex_malloc,
            "Not enough memory for chunkbuf");

    big_array_init(&chunk_array, chunkbuf, bb->dtype, 2, dims, NULL);
    big_array_iter_init(&array_iter, array);

    toread = array->size / bb->nmemb;

    ptrdiff_t abs = bb->foffset[ptr->fileid] + ptr->roffset + toread;
    RAISEIF(abs > bb->size,
                ex_eof,
                "Reading beyond the block `%s` at (%d:%td)",
                bb->basename, ptr->fileid, ptr->roffset * felsize);

    while(toread > 0 && ! big_block_eof(bb, ptr)) {
        size_t chunk_size = CHUNK_SIZE;
        /* remaining items in the file */
        if(chunk_size > bb->fsize[ptr->fileid] - ptr->roffset) {
            chunk_size = bb->fsize[ptr->fileid] - ptr->roffset;
        }
        /* remaining items to read */
        if(chunk_size > toread) {
            chunk_size = toread;
        }
        RAISEIF(chunk_size == 0,
            ex_insuf,
            "Insufficient number of items in file `%s' at (%d:%td)",
            bb->basename, ptr->fileid, ptr->roffset * felsize);

        /* read to the beginning of chunk */
        big_array_iter_init(&chunk_iter, &chunk_array);

        fp = _big_file_open_a_file(bb->basename, ptr->fileid, "r");
        RAISEIF(fp == NULL,
                ex_open,
                NULL);
        RAISEIF(0 > fseek(fp, ptr->roffset * felsize, SEEK_SET),
                ex_seek,
                "Failed to seek in block `%s' at (%d:%td) (%s)", 
                bb->basename, ptr->fileid, ptr->roffset * felsize, strerror(errno));
        RAISEIF(chunk_size != fread(chunkbuf, felsize, chunk_size, fp),
                ex_read,
                "Failed to read in block `%s' at (%d:%td) (%s)",
                bb->basename, ptr->fileid, ptr->roffset * felsize, strerror(errno));
        fclose(fp);

        /* now translate the data from chunkbuf to mptr */
        _dtype_convert(&array_iter, &chunk_iter, chunk_size * bb->nmemb);

        toread -= chunk_size;
        RAISEIF(0 != big_block_seek_rel(bb, ptr, chunk_size),
                ex_blockseek,
                NULL);
    }

    free(chunkbuf);
    return 0;
ex_insuf:
ex_read:
ex_seek:
    fclose(fp);
ex_blockseek:
ex_open:
ex_eof:
    free(chunkbuf);
ex_malloc:
    return -1;
}

int
big_block_write(BigBlock * bb, BigBlockPtr * ptr, BigArray * array)
{
    if(array->size == 0) return 0;
    /* the file header is modified */
    bb->dirty = 1;
    char * chunkbuf = malloc(CHUNK_BYTES);
    int felsize = dtype_itemsize(bb->dtype) * bb->nmemb;
    size_t CHUNK_SIZE = CHUNK_BYTES / felsize;

    BigArray chunk_array = {0};
    size_t dims[2];
    dims[0] = CHUNK_SIZE;
    dims[1] = bb->nmemb;

    BigArrayIter chunk_iter;
    BigArrayIter array_iter;
    ptrdiff_t towrite = 0;
    FILE * fp;

    RAISEIF(chunkbuf == NULL,
            ex_malloc,
            "not enough memory for chunkbuf of size %d bytes", CHUNK_BYTES);

    big_array_init(&chunk_array, chunkbuf, bb->dtype, 2, dims, NULL);
    big_array_iter_init(&array_iter, array);

    towrite = array->size / bb->nmemb;

    ptrdiff_t abs = bb->foffset[ptr->fileid] + ptr->roffset + towrite;
    RAISEIF(abs > bb->size,
                ex_eof,
                "Writing beyond the block `%s` at (%d:%td)",
                bb->basename, ptr->fileid, ptr->roffset * felsize);

    while(towrite > 0 && ! big_block_eof(bb, ptr)) {
        size_t chunk_size = CHUNK_SIZE;
        /* remaining items in the file */
        if(chunk_size > bb->fsize[ptr->fileid] - ptr->roffset) {
            chunk_size = bb->fsize[ptr->fileid] - ptr->roffset;
        }
        /* remaining items to read */
        if(chunk_size > towrite) {
            chunk_size = towrite;
        }
        /* write from the beginning of chunk */
        big_array_iter_init(&chunk_iter, &chunk_array);

        /* now translate the data to format in the file*/
        _dtype_convert(&chunk_iter, &array_iter, chunk_size * bb->nmemb);

        sysvsum(&bb->fchecksum[ptr->fileid], chunkbuf, chunk_size * felsize);

        fp = _big_file_open_a_file(bb->basename, ptr->fileid, "r+");
        RAISEIF(fp == NULL,
                ex_open,
                NULL);
        RAISEIF(0 > fseek(fp, ptr->roffset * felsize, SEEK_SET),
                ex_seek,
                "Failed to seek in block `%s' at (%d:%td) (%s)", 
                bb->basename, ptr->fileid, ptr->roffset * felsize, strerror(errno));
        RAISEIF(chunk_size != fwrite(chunkbuf, felsize, chunk_size, fp),
                ex_write,
                "Failed to write in block `%s' at (%d:%td) (%s)",
                bb->basename, ptr->fileid, ptr->roffset * felsize, strerror(errno));
        fclose(fp);

        towrite -= chunk_size;
        RAISEIF(0 != big_block_seek_rel(bb, ptr, chunk_size),
                ex_blockseek, NULL);
    }
    free(chunkbuf);
    return 0;
ex_write:
ex_seek:
    fclose(fp);
ex_open:
ex_blockseek:
ex_eof:
    free(chunkbuf);
ex_malloc:
    return -1;
}

/**
 * dtype stuff
 * */

#define MACHINE_ENDIANNESS MACHINE_ENDIAN_F()
static char
MACHINE_ENDIAN_F(void)
{
    uint32_t i = 0x01234567;
    return ((*((uint8_t*)(&i))) == 0x67)?'<':'>';
}

int
dtype_normalize(char * dst, const char * src)
{
/* normalize a dtype, so that
 * dst[0] is the endian-ness
 * dst[1] is the type kind char
 * dtype[2:] is the width
 * */
    memset(dst, 0, 8);
    switch(src[0]) {
        case '<':
        case '>':
        case '=':
            strcpy(dst, src);
        break;
        default:
            dst[0] = '=';
            strcpy(dst + 1, src);
    }
    if(dst[0] == '=') {
        dst[0] = MACHINE_ENDIANNESS;
    }
    return 0;
}

int
dtype_itemsize(const char * dtype)
{
    char ndtype[8];
    dtype_normalize(ndtype, dtype);
    return atoi(&ndtype[2]);
}

int
big_array_init(BigArray * array, void * buf, const char * dtype, int ndim, const size_t dims[], const ptrdiff_t strides[])
{

    memset(array, 0, sizeof(array[0]));

    dtype_normalize(array->dtype, dtype);
    array->data = buf;
    array->ndim = ndim;
    int i;
    memset(array->dims, 0, sizeof(ptrdiff_t) * 32);
    memset(array->strides, 0, sizeof(ptrdiff_t) * 32);
    array->size = 1;
    for(i = 0; i < ndim; i ++) {
        array->dims[i] = dims[i];
        array->size *= dims[i];
    }
    if(strides != NULL) {
        for(i = 0; i < ndim; i ++) {
            array->strides[i] = strides[i];
        }
    } else {
        array->strides[ndim - 1] = dtype_itemsize(dtype);
        for(i = ndim - 2; i >= 0; i --) {
            array->strides[i] = array->strides[i + 1] * array->dims[i + 1];
        }
    }
    return 0;
}

int
big_array_iter_init(BigArrayIter * iter, BigArray * array)
{
    memset(iter, 0, sizeof(iter[0]));

    iter->array = array;

    memset(iter->pos, 0, sizeof(ptrdiff_t) * 32);
    iter->dataptr = array->data;

    /* see if the iter is contiguous */
    size_t elsize = dtype_itemsize(array->dtype);

    int i = 0; 
    ptrdiff_t stride_contiguous = elsize;
    iter->contiguous = 1;
    for(i = array->ndim - 1; i >= 0; i --) {
        if(array->strides[i] != stride_contiguous) {
            iter->contiguous = 0;
            break;
        }
        stride_contiguous *= array->dims[i];
    }
    return 0;
}

void
big_array_iter_advance(BigArrayIter * iter)
{
    BigArray * array = iter->array;

    if(iter->contiguous) {
        iter->dataptr = (char*) iter->dataptr + array->strides[array->ndim - 1];
        return;
    }
    int k;
    iter->pos[array->ndim - 1] ++;
    iter->dataptr = ((char*) iter->dataptr) + array->strides[array->ndim - 1];
    for(k = array->ndim - 1; k >= 0; k --) {
        if(iter->pos[k] == array->dims[k]) {
            iter->dataptr = ((char*) iter->dataptr) - array->strides[k] * iter->pos[k];
            iter->pos[k] = 0;
            if(k > 0) {
                iter->pos[k - 1] ++;
                iter->dataptr = ((char*) iter->dataptr) + array->strides[k - 1];
            }
        } else {
            break;
        }
    }
}

/* format data in dtype to a string in buffer */
void
dtype_format(char * buffer, const char * dtype, const void * data, const char * fmt)
{
    char ndtype[8];
    char ndtype2[8];
    union {
        char *S1;
        int64_t *i8;
        uint64_t *u8;
        double *f8;
        int32_t *i4;
        uint32_t *u4;
        float *f4;
        void * v;
    } p;

    /* handle the endianness stuff in case it is not machine */
    char converted[128];

    dtype_normalize(ndtype2, dtype);
    ndtype2[0] = '=';
    dtype_normalize(ndtype, ndtype2);
    dtype_convert_simple(converted, ndtype, data, dtype, 1);

    p.v = converted;
#define _HANDLE_FMT_(dtype, defaultfmt) \
    if(!strcmp(ndtype + 1, # dtype)) { \
        if(fmt == NULL) fmt = defaultfmt; \
        sprintf(buffer, fmt, *p.dtype); \
    } else

    _HANDLE_FMT_(S1, "%c")
    _HANDLE_FMT_(i8, "%ld")
    _HANDLE_FMT_(i4, "%d")
    _HANDLE_FMT_(u8, "%lu")
    _HANDLE_FMT_(u4, "%u")
    _HANDLE_FMT_(f8, "%g")
    _HANDLE_FMT_(f4, "%g")
    abort();
#undef _HANDLE_FMT_
}

/* parse data in dtype to a string in buffer */
void
dtype_parse(const char * buffer, const char * dtype, void * data, const char * fmt)
{
    char ndtype[8];
    char ndtype2[8];
    union {
        char *S1;
        int64_t *i8;
        uint64_t *u8;
        double *f8;
        int32_t *i4;
        uint32_t *u4;
        float *f4;
        void * v;
    } p;

    /* handle the endianness stuff in case it is not machine */
    char converted[128];

    dtype_normalize(ndtype2, dtype);
    ndtype2[0] = '=';
    dtype_normalize(ndtype, ndtype2);

    p.v = converted;
#define _HANDLE_FMT_(dtype, defaultfmt) \
    if(!strcmp(ndtype + 1, # dtype)) { \
        if(fmt == NULL) fmt = defaultfmt; \
        sscanf(buffer, fmt, p.dtype); \
    } else

    _HANDLE_FMT_(S1, "%c")
    _HANDLE_FMT_(i8, "%ld")
    _HANDLE_FMT_(i4, "%d")
    _HANDLE_FMT_(u8, "%lu")
    _HANDLE_FMT_(u4, "%u")
    _HANDLE_FMT_(f8, "%lf")
    _HANDLE_FMT_(f4, "%f")
    abort();
#undef _HANDLE_FMT_

    dtype_convert_simple(data, dtype, converted, ndtype, 1);

}

static int
dtype_convert_simple(void * dst, const char * dstdtype, const void * src, const char * srcdtype, size_t nmemb)
{
    BigArray dst_array, src_array;
    BigArrayIter dst_iter, src_iter;
    big_array_init(&dst_array, dst, dstdtype, 1, &nmemb, NULL);
    big_array_init(&src_array, (void*) src, srcdtype, 1, &nmemb, NULL);
    big_array_iter_init(&dst_iter, &dst_array);
    big_array_iter_init(&src_iter, &src_array);
    return _dtype_convert(&dst_iter, &src_iter, nmemb);
}

static void cast(BigArrayIter * dst, BigArrayIter * src, size_t nmemb);
static void byte_swap(BigArrayIter * array, size_t nmemb);
int
_dtype_convert(BigArrayIter * dst, BigArrayIter * src, size_t nmemb)
{
    /* cast buf2 of dtype2 into buf1 of dtype1 */
    /* match src to machine endianness */
    if(src->array->dtype[0] != MACHINE_ENDIANNESS) {
        BigArrayIter iter = *src;
        byte_swap(&iter, nmemb);
    }

    BigArrayIter iter1 = *dst;
    BigArrayIter iter2 = *src;
    cast(&iter1, &iter2, nmemb);

    /* match dst to machine endianness */
    if(dst->array->dtype[0] != MACHINE_ENDIANNESS) {
        BigArrayIter iter = *dst;
        byte_swap(&iter, nmemb);
    }
    *dst = iter1;
    *src = iter2;
    return 0;
}


static void
byte_swap(BigArrayIter * iter, size_t nmemb)
{
    /* swap a buffer in-place */
    int elsize = dtype_itemsize(iter->array->dtype);
    if(elsize == 1) return;
    /* need byte swap; do it now on buf2 */
    ptrdiff_t i;
    int half = elsize << 1;
    for(i = 0; i < nmemb; i ++) {
        int j;
        char * ptr = iter->dataptr;
        for(j = 0; j < half; j ++) {
            char tmp = ptr[j];
            ptr[j] = ptr[elsize - j - 1];
            ptr[elsize - j - 1] = tmp;
        }
        big_array_iter_advance(iter);
    }
}

#define CAST_CONVERTER(d1, t1, d2, t2) \
if(!strcmp(d1, dst->array->dtype + 1) && !strcmp(d2, src->array->dtype + 1)) { \
    ptrdiff_t i; \
    for(i = 0; i < nmemb; i ++) { \
        t1 * p1 = dst->dataptr; t2 * p2 = src->dataptr; \
        * p1 = * p2; \
        big_array_iter_advance(dst); big_array_iter_advance(src); \
    } \
    return; \
}
static void
cast(BigArrayIter * dst, BigArrayIter * src, size_t nmemb)
{
    /* doing cast assuming native byte order */

    /* convert buf2 to buf1, both are native;
     * dtype has no endian-ness prefix
     *   */
    if((dst->contiguous && src->contiguous) 
    && !strcmp(dst->array->dtype + 1, src->array->dtype + 1)) {
        memcpy(dst->dataptr, src->dataptr, nmemb * dst->array->strides[dst->array->ndim-1]);
        dst->dataptr = (char*) dst->dataptr + nmemb * dst->array->strides[dst->array->ndim - 1];
        src->dataptr = (char*) src->dataptr + nmemb * src->array->strides[src->array->ndim - 1];
        return;
    }
    if(!strcmp(dst->array->dtype + 1, "i8")) {
        CAST_CONVERTER("i8", int64_t, "i8", int64_t);
        CAST_CONVERTER("i8", int64_t, "i4", int32_t);
        CAST_CONVERTER("i8", int64_t, "u4", uint32_t);
        CAST_CONVERTER("i8", int64_t, "u8", uint64_t);
        CAST_CONVERTER("i8", int64_t, "f8", double);
        CAST_CONVERTER("i8", int64_t, "f4", float);
    } else
    if(!strcmp(dst->array->dtype + 1, "u8")) {
        CAST_CONVERTER("u8", uint64_t, "u8", uint64_t);
        CAST_CONVERTER("u8", uint64_t, "u4", uint32_t);
        CAST_CONVERTER("u8", uint64_t, "i4", int32_t);
        CAST_CONVERTER("u8", uint64_t, "i8", int64_t);
        CAST_CONVERTER("u8", uint64_t, "f8", double);
        CAST_CONVERTER("u8", uint64_t, "f4", float);
    } else 
    if(!strcmp(dst->array->dtype + 1, "f8")) {
        CAST_CONVERTER("f8", double, "f8", double);
        CAST_CONVERTER("f8", double, "f4", float);
        CAST_CONVERTER("f8", double, "i4", int32_t);
        CAST_CONVERTER("f8", double, "i8", int64_t);
        CAST_CONVERTER("f8", double, "u4", uint32_t);
        CAST_CONVERTER("f8", double, "u8", uint64_t);
    } else
    if(!strcmp(dst->array->dtype + 1, "i4")) {
        CAST_CONVERTER("i4", int32_t, "i4", int32_t);
        CAST_CONVERTER("i4", int32_t, "i8", int64_t);
        CAST_CONVERTER("i4", int32_t, "u4", uint32_t);
        CAST_CONVERTER("i4", int32_t, "u8", uint64_t);
        CAST_CONVERTER("i4", int32_t, "f8", double);
        CAST_CONVERTER("i4", int32_t, "f4", float);
    } else
    if(!strcmp(dst->array->dtype + 1, "u4")) {
        CAST_CONVERTER("u4", uint32_t, "u4", uint32_t);
        CAST_CONVERTER("u4", uint32_t, "u8", uint64_t);
        CAST_CONVERTER("u4", uint32_t, "i4", int32_t);
        CAST_CONVERTER("u4", uint32_t, "i8", int64_t);
        CAST_CONVERTER("u4", uint32_t, "f8", double);
        CAST_CONVERTER("u4", uint32_t, "f4", float);
    } else
    if(!strcmp(dst->array->dtype + 1, "f4")) {
        CAST_CONVERTER("f4", float, "f4", float);
        CAST_CONVERTER("f4", float, "f8", double);
        CAST_CONVERTER("f4", float, "i4", int32_t);
        CAST_CONVERTER("f4", float, "i8", int64_t);
        CAST_CONVERTER("f4", float, "u4", uint32_t);
        CAST_CONVERTER("f4", float, "u8", uint64_t);
    } else
    if(!strcmp(dst->array->dtype + 1, "S1")) {
        CAST_CONVERTER("S1", char, "S1", char);
    }
    /* */
    fprintf(stderr, "Unsupported conversion from %s to %s\n", src->array->dtype, dst->array->dtype);
    abort();
}
#undef CAST_CONVERTER

static void
sysvsum(unsigned int * sum, void * buf, size_t size)
{
    unsigned int thisrun = *sum;
    unsigned char * cp = buf;
    while(size --)    
        thisrun += *(cp++);
    *sum = thisrun;
}

/*
 * Internal API for AttrSet objects;
 * */

static int
attrset_read_attr_set_v1(BigAttrSet * attrset, const char * basename)
{
    attrset->dirty = 0;

    FILE * fattr = _big_file_open_a_file(basename, FILEID_ATTR, "r");
    if(fattr == NULL) {
        big_file_clear_error_message();
        return 0;
    }
    int nmemb;
    int lname;
    char dtype[8];
    char * data;
    char * name;
    while(!feof(fattr)) {
        if(1 != fread(&nmemb, sizeof(int), 1, fattr)) break;
        RAISEIF(
            (1 != fread(&lname, sizeof(int), 1, fattr)) ||
            (1 != fread(&dtype, 8, 1, fattr)),
            ex_fread,
            "Failed to read from file"
                )
        int ldata = dtype_itemsize(dtype) * nmemb;
        data = alloca(ldata);
        name = alloca(lname + 1);
        RAISEIF(
            (1 != fread(name, lname, 1, fattr)) ||
            (1 != fread(data, ldata, 1, fattr)),
            ex_fread,
            "Failed to read from file");

        name[lname] = 0;
        RAISEIF(0 != attrset_set_attr(attrset, name, data, dtype, nmemb),
            ex_set_attr,
            NULL);
    } 
    attrset->dirty = 0;
    fclose(fattr);
    return 0;
ex_set_attr:
ex_fread:
    attrset->dirty = 0;
    fclose(fattr);
    return -1;
}

static int _isblank(int ch) {
    return ch == ' ' || ch == '\t';
}

static int
attrset_read_attr_set_v2(BigAttrSet * attrset, const char * basename)
{
    attrset->dirty = 0;

    FILE * fattr = _big_file_open_a_file(basename, FILEID_ATTR_V2, "r");
    if(fattr == NULL) {
        big_file_clear_error_message();
        return 0;
    }
    fseek(fattr, 0, SEEK_END);
    long size = ftell(fattr);
    char * buffer = (char*) malloc(size + 1);
    unsigned char * data = (unsigned char * ) malloc(size + 1);
    fseek(fattr, 0, SEEK_SET);
    RAISEIF(size != fread(buffer, 1, size, fattr),
            ex_read_file,
            "Failed to read attribute file\n");
    fclose(fattr);
    buffer[size] = 0;

    /* now parse the v2 attr file.*/
    long i = 0;
    #define ATTRV2_EXPECT(variable) while(_isblank(buffer[i])) i++; \
                    char * variable = buffer + i; \
                    while(!_isblank(buffer[i])) i++; buffer[i] = 0; i++;
    while(buffer[i]) {
        ATTRV2_EXPECT(name);
        ATTRV2_EXPECT(dtype);
        ATTRV2_EXPECT(rawlength);
        ATTRV2_EXPECT(rawdata);
        /* skip the reset of the line */
        while(buffer[i] != '\n' && buffer[i]) i ++;
        if(buffer[i] == '\n') i++;

        int nmemb = atoi(rawlength);
        int itemsize = dtype_itemsize(dtype);

        RAISEIF(nmemb * itemsize * 2!= strlen(rawdata),
            ex_parse_attr,
            "NMEMB and data mismiatch: %d x %d (%s) * 2 != %d",
            nmemb, itemsize, dtype, strlen(rawdata));

        int j, k;
        for(k = 0, j = 0; k < nmemb * itemsize; k ++, j += 2) {
            char buf[3];
            buf[0] = rawdata[j];
            buf[1] = rawdata[j+1];
            buf[2] = 0;
            unsigned int byte = strtoll(buf, NULL, 16);
            data[k] = byte;
        }
        RAISEIF(0 != attrset_set_attr(attrset, name, data, dtype, nmemb),
            ex_set_attr,
            NULL);
    } 
    free(data);
    free(buffer);
    attrset->dirty = 0;
    return 0;

ex_read_file:
ex_parse_attr:
ex_set_attr:
    attrset->dirty = 0;
    free(buffer);
    free(data);
    return -1;
}
static int
attrset_write_attr_set_v2(BigAttrSet * attrset, const char * basename)
{
    static char conv[] = "0123456789ABCDEF";
    attrset->dirty = 0;

    FILE * fattr = _big_file_open_a_file(basename, FILEID_ATTR_V2, "w");
    RAISEIF(fattr == NULL,
            ex_open,
            NULL);

    ptrdiff_t i;
    for(i = 0; i < attrset->listused; i ++) {
        BigAttr * a = & attrset->attrlist[i];
        int itemsize = dtype_itemsize(a->dtype);
        int ldata = itemsize * a->nmemb;

        char * rawdata = malloc(ldata * 2 + 1);
        char * textual = malloc(a->nmemb * 32 + 1);
        textual[0] = 0;
        unsigned char * adata = (unsigned char*) a->data;
        int j, k; 
        for(j = 0, k = 0; k < ldata; k ++, j+=2) {
            rawdata[j] = conv[adata[k] / 16];
            rawdata[j + 1] = conv[adata[k] % 16];
        }
        rawdata[j] = 0;
        for(j = 0; j < a->nmemb; j ++) {
            if(a->dtype[1] != 'S') {
                char buf[128];
                dtype_format(buf, a->dtype, &adata[j * itemsize], NULL);
                strcat(textual, buf);
                if(j != a->nmemb - 1) {
                    strcat(textual, " ");
                }
            } else {
                char buf[] = {adata[j], 0};
                if(buf[0] == '\n') {
                    strcat(textual, "...");
                    break;
                } if(buf[0] == 0) {
                    break;
                }
                strcat(textual, buf);
            }
        }
        int rt = fprintf(fattr, "%s %s %d %s #HUMANE [ %s ]\n", 
                a->name, a->dtype, a->nmemb, rawdata, textual
               );
        free(rawdata);
        free(textual);
        RAISEIF(rt <= 0,
            ex_write,
            "Failed to write to file");
    } 
    fclose(fattr);
    return 0;
ex_write:
    fclose(fattr);
ex_open:
    return -1;
}

static int attr_cmp(const void * p1, const void * p2) {
    const BigAttr * c1 = p1;
    const BigAttr * c2 = p2;
    return strcmp(c1->name, c2->name);
}

static BigAttr *
attrset_append_attr(BigAttrSet * attrset)
{
    while(attrset->listsize - attrset->listused < 1) {
        attrset->attrlist = realloc(attrset->attrlist, attrset->listsize * 2 * sizeof(BigAttr));
        attrset->listsize *= 2;
    }
    BigAttr * a = & (attrset->attrlist[attrset->listused++]);
    memset(a, 0, sizeof(BigAttr));
    return a;
}

static int
attrset_add_attr(BigAttrSet * attrset, const char * attrname, const char * dtype, int nmemb)
{
    size_t size = dtype_itemsize(dtype) * nmemb + strlen(attrname) + 1;
    while(attrset->bufsize - attrset->bufused < size) {
        int i;
        for(i = 0; i < attrset->listused; i ++) {
            attrset->attrlist[i].data -= (ptrdiff_t) attrset->attrbuf;
            attrset->attrlist[i].name -= (ptrdiff_t) attrset->attrbuf;
        }
        attrset->attrbuf = realloc(attrset->attrbuf, attrset->bufsize * 2);
        attrset->bufsize *= 2;
        for(i = 0; i < attrset->listused; i ++) {
            attrset->attrlist[i].data += (ptrdiff_t) attrset->attrbuf;
            attrset->attrlist[i].name += (ptrdiff_t) attrset->attrbuf;
        }
    }
    char * free = attrset->attrbuf + attrset->bufused;
    attrset->bufused += size;

    BigAttr * n = attrset_append_attr(attrset);

    n->nmemb = nmemb;
    memset(n->dtype, 0, 8);
    dtype_normalize(n->dtype, dtype);

    n->name = free;
    strcpy(free, attrname);
    free += strlen(attrname) + 1;
    n->data = free;

    qsort(attrset->attrlist, attrset->listused, sizeof(BigAttr), attr_cmp);
    return 0;
}

static BigAttr *
attrset_lookup_attr(BigAttrSet * attrset, const char * attrname)
{
    BigAttr lookup = {0};
    lookup.name = (char*) attrname;
    BigAttr * found = bsearch(&lookup, attrset->attrlist, attrset->listused, sizeof(BigAttr), attr_cmp);
    return found;
}

static int
attrset_remove_attr(BigAttrSet * attrset, const char * attrname)
{
    BigAttr *attr = attrset_lookup_attr(attrset, attrname);
    if(attr) {
        ptrdiff_t ind = attr - attrset->attrlist;
        memmove(&attrset->attrlist[ind], &attrset->attrlist[ind + 1],
            (attrset->listused - ind - 1) * sizeof(BigAttr));
        attrset->listused -= 1;
    }
    return 0;
}

static BigAttr *
attrset_list_attrs(BigAttrSet * attrset, size_t * count)
{
    *count = attrset->listused;
    return attrset->attrlist;
}

static int
attrset_set_attr(BigAttrSet * attrset, const char * attrname, const void * data, const char * dtype, int nmemb)
{
    BigAttr * attr;
    attrset->dirty = 1;

    attrset_remove_attr(attrset, attrname);
    /* add ensures the dtype has been normalized! */
    RAISEIF(0 != attrset_add_attr(attrset, attrname, dtype, nmemb),
            ex_add,
            "Failed to add attr");
    attr = attrset_lookup_attr(attrset, attrname);
    RAISEIF(attr->nmemb != nmemb,
            ex_mismatch,
            "attr nmemb mismatch");
    dtype_convert_simple(attr->data, attr->dtype, data, dtype, attr->nmemb);
    return 0;

ex_mismatch:
ex_add:
    return -1;
}

static int
attrset_get_attr(BigAttrSet * attrset, const char * attrname, void * data, const char * dtype, int nmemb)
{
    BigAttr * found = attrset_lookup_attr(attrset, attrname);
    RAISEIF(!found, ex_notfound, "attr not found");
    RAISEIF(found->nmemb != nmemb, ex_mismatch, "attr nmemb mismatch");
    dtype_convert_simple(data, dtype, found->data, found->dtype, found->nmemb);
    return 0;

ex_mismatch:
ex_notfound:
    return -1;
}

static BigAttrSet *
attrset_create()
{
    BigAttrSet * attrset = calloc(1, sizeof(BigAttrSet));
    attrset->attrbuf = malloc(128);
    attrset->bufsize = 128;
    attrset->bufused = 0;
    attrset->attrlist = malloc(sizeof(BigAttr) * 16);
    attrset->listsize = 16;
    attrset->listused = 0;

    return attrset;
}

static void
attrset_free(BigAttrSet * attrset)
{
    free(attrset->attrbuf);
    free(attrset->attrlist);
    free(attrset);
}

void big_attrset_set_dirty(BigAttrSet * attrset, int dirty)
{
    attrset->dirty = dirty;
}
void * big_attrset_pack(BigAttrSet * attrset, size_t * bytes)
{
    size_t n = 0;
    n += sizeof(BigAttrSet);
    n += attrset->bufused;
    n += sizeof(BigAttr) * attrset->listused;
    char * buf = calloc(n, 1);
    char * p = buf;
    char * attrbuf = (char*) (p + sizeof(BigAttrSet));
    BigAttr * attrlist = (BigAttr *) (attrbuf + attrset->bufused);
    memcpy(p, attrset, sizeof(BigAttrSet));
    memcpy(attrbuf, attrset->attrbuf, attrset->bufused);
    memcpy(attrlist, attrset->attrlist, attrset->listused * sizeof(BigAttr));
    int i = 0;
    for(i = 0; i < attrset->listused; i ++) {
        attrlist[i].data -= (ptrdiff_t) attrset->attrbuf;
        attrlist[i].name -= (ptrdiff_t) attrset->attrbuf;
    }

    *bytes = n;

    return (void*) p;
}

BigAttrSet * big_attrset_unpack(void * p)
{
    BigAttrSet * attrset = calloc(1, sizeof(attrset[0]));
    memcpy(attrset, p, sizeof(BigAttrSet));
    p += sizeof(BigAttrSet);
    attrset->attrbuf = malloc(attrset->bufsize);
    attrset->attrlist = malloc(attrset->listsize * sizeof(BigAttr));
    memcpy(attrset->attrbuf, p, attrset->bufused);
    p += attrset->bufused;
    memcpy(attrset->attrlist, p, attrset->listused * sizeof(BigAttr));
    int i = 0;
    for(i = 0; i < attrset->listused; i ++) {
        attrset->attrlist[i].data += (ptrdiff_t) attrset->attrbuf;
        attrset->attrlist[i].name += (ptrdiff_t) attrset->attrbuf;
    }
    return attrset;
}

/* File Path */

FILE *
_big_file_open_a_file(const char * basename, int fileid, char * mode)
{
    char * filename;
    if(fileid == FILEID_HEADER) {
        filename = _path_join(basename, EXT_HEADER);
    } else
    if(fileid == FILEID_ATTR) {
        filename = _path_join(basename, EXT_ATTR);
    } else
    if(fileid == FILEID_ATTR_V2) {
        filename = _path_join(basename, EXT_ATTR_V2);
    } else {
        char d[128];
        sprintf(d, EXT_DATA, fileid);
        filename = _path_join(basename, d);
    }
    FILE * fp = fopen(filename, mode);
    RAISEIF(fp == NULL,
        ex_fopen,
        "Failed to open physical file `%s' with mode `%s' (%s)",
        filename, mode, strerror(errno));
ex_fopen:
    free(filename);
    return fp;
}
/* make subdir rel to pathname, recursively making parents */
int
_big_file_mksubdir_r(const char * pathname, const char * subdir)
{
    char * subdirname = strdup(subdir);
    char * p;

    char * mydirname;
    for(p = subdirname; *p; p ++) {
        if(*p != '/') continue;
        *p = 0;
        mydirname = _path_join(pathname, subdirname);
        mkdir(mydirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        free(mydirname);
        *p = '/';
    }
    mydirname = _path_join(pathname, subdirname);
    mkdir(mydirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    struct stat buf;
    RAISEIF(0 != stat(mydirname, &buf),
            ex_mkdir,
            "Failed to create directory structure at `%s' (%s)", 
            mydirname,
            strerror(errno)
    );
    free(subdirname);
    free(mydirname);
    return 0;
ex_mkdir:
    free(mydirname);
    return -1;
}

