#ifndef  __HEAP_H__
#define __HEAP_H___
int heap_init(size_t allocation);
void * heap_allocate(size_t bytes);
int heap_return0(void * ptr);
#define heap_return(x) (heap_return0(x), x = NULL)
size_t heap_get_max_usage();
size_t heap_get_free_bytes();
size_t heap_get_total_bytes();
#endif
