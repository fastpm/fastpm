#ifndef TIMER_H
#define TIMER_H 1

enum Category {INIT, LPT, STEPPING, SNP};

void timer_set_category(enum Category new_cat);
void timer_start(const char * sub);
void timer_stop(const char * sub);
void timer_print();

#endif
