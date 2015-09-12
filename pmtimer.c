/****
 * Profiling library
 * 
 * from Jun Koda's PM framework utilities.
 *****/
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "msg.h"
#include "pmtimer.h"

#define nCategory 4
#define nSubCategory 15

static const char * CatName[]= {"Init", "2LPT", "Stepping", "Snapshot"};

static int initialized= 0;
static enum Category Cat;
static double Time[nCategory][nSubCategory], tBegin[nCategory][nSubCategory];
static char * Name[nCategory][nSubCategory];

static int timer_find(const char * sub) {
    int isub;
    for(isub = 0; isub < nSubCategory; isub ++) {
        if(Name[Cat][isub] == NULL) {
            Name[Cat][isub] = strdup(sub);
            break;
        }
        if(!strcmp(Name[Cat][isub], sub)) break;
    }
    if(isub == nSubCategory) {
        msg_abort(-1, "Too many sub timers. Increase nSubCategory in pmtimer.c and recompile.");
    }
    return isub;
    
}

static double now()
{
    struct timeval tp;
    //struct timezone tzp;
    //gettimeofday(&tp,&tzp);
    gettimeofday(&tp, 0);

    return (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6;
}

void timer_init()
{
    int i, j;
    for(i=0; i<nCategory; i++) {
        for(j=0; j<nSubCategory; j++) {
            Time[i][j]= 0.0;
            tBegin[i][j]= 0.0;
            Name[i][j] = NULL;
        }
        Cat = i;
        timer_find("total");
    }

    initialized= 1;
}

void timer_set_category(enum Category new_cat)
{
    double nw= now();

    if(initialized)
        Time[Cat][0] += nw - tBegin[Cat][0];
    else 
        timer_init();

    Cat= new_cat;
    tBegin[Cat][0]= nw;
}

void timer_start(const char * sub)
{
    int isub = timer_find(sub);
    tBegin[Cat][isub]= now();
}

void timer_stop(const char * sub)
{
    int isub = timer_find(sub);
    Time[Cat][isub] += now() - tBegin[Cat][isub];
}

void timer_print()
{
    timer_set_category(0);
    double total= 0.0;
    int i, icat, isub;
    for(i=0; i<nCategory; i++)
        total += Time[i][0];

    for(icat=0; icat<nCategory; icat++) {
        msg_printf(info, "%-16s %7.2f   %4.1f%%\n", 
                CatName[icat], Time[icat][0], 100.0*Time[icat][0]/total);

        for(isub=1; isub<nSubCategory; isub++) {
            if(Time[icat][isub] > 0.0 && Name[icat][isub] != NULL )
                msg_printf(info, "  %-14s %7.2f   %4.1f%%\n", 
                        Name[icat][isub], Time[icat][isub], 100*Time[icat][isub]/total);
        }
    }
    msg_printf(info, "----------------------------------\n");
    msg_printf(info, "%-16s %7.2f sec\n", "Total", total);

}


