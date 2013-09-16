#ifndef TIMER_H
#define TIMER_H 1

enum Category {Init, LPT, COLA, Snp};
enum SubCategory {all, fft, assign, force_mesh, pforce, check, comm, evolve, write, kd_build, kd_link, interp, global, sub};

void timer_set_category(enum Category new_cat);
void timer_start(enum SubCategory sub);
void timer_stop(enum SubCategory sub);
void timer_print();

#endif
