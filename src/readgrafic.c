#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

void 
read_grafic_gaussian(PM * pm, FastPMFloat * g_x, char * filename)
{
    ptrdiff_t ind;
    int d;
    ptrdiff_t i[3] = {0, 0, 0};

    FILE * fp = fopen(filename, "r");

    struct {
        int32_t bs1;
        int32_t n[3];
        int32_t seed;
        int32_t bs2;
    } header;

    if(1 != fread(&header, sizeof(header), 1, fp)) {
        fastpm_raise(-1, "file not in BigMD noise format.\n");
    }

    ind = 0;

    if(header.bs1 != 16) 
        fastpm_raise(-1, "file not in BigMD noise format.\n");

    for(d = 0; d < 3; d ++) {
        /* BigMD is in */
        int permute[3] = {2, 1, 0};
        if(header.n[d] != pm_nmesh(pm)[permute[d]]) {
            fastpm_raise(-1, "file is in %d, but simulation is in %d.\n",
                header.n[d], pm_nmesh(pm)[permute[d]]);
        }
    }

    float * buf = malloc(sizeof(float) * header.n[0] * pm_i_region(pm)->size[1]);

    for(i[0] = 0; i[0] < pm_i_region(pm)->size[0]; i[0] ++) {
        ptrdiff_t i_abs[3];

        i_abs[0] = i[0] + pm_i_region(pm)->start[0];

        ptrdiff_t offset = (/*header*/16 + 8) + 
                    i_abs[0] * (header.n[0] * header.n[1] * 4 + 8); 
        int32_t bs;

        fseek(fp, offset, SEEK_SET);

        if(1 != fread(&bs, 4, 1, fp))
            fastpm_raise(-1, "file not in BigMD noise format.\n");

        if(bs != 4 * header.n[0] * header.n[1])
            fastpm_raise(-1, "file size is wrong\n");

        fseek(fp, pm_i_region(pm)->start[1] * header.n[1] * 4, SEEK_CUR);

        if(header.n[0] * pm_i_region(pm)->size[1] != fread(buf, 4, header.n[0] * pm_i_region(pm)->size[1], fp))
            fastpm_raise(-1, "file not in BigMD noise format.\n");

        int p = 0;
        for(i[1] = 0; i[1] < pm_i_region(pm)->size[1]; i[1] ++) {
            /* note that size[2] is Nmesh[2] == n[0] and start[2] is 0 */
            for(i[2] = 0; i[2] < pm_i_region(pm)->size[2]; i[2] ++) {
                ind = 0;
                for(d = 0; d < 3; d++) {
                    ind += i[d] * pm_i_region(pm)->strides[d];
                }
                g_x[ind] = buf[p];
                p ++;
            }
        }
    }
    free(buf);
}

