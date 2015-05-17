#ifndef  __HEAP_H__
#define __HEAP_H___
int heap_init(size_t allocation);
void * heap_allocate(size_t bytes);
int heap_return(void * ptr);
size_t heap_get_max_usage();
#endif
