typedef double (*kd_gravity_func)( double r, void * userdata);
double kd_grav(KDNode * node, double * pos, double openangle,
        kd_gravity_func gravity, 
        void * userdata) {
    KDStore * t0 = node->store;
    int Nd = node->store->input.dims[1];
    int d;
    double *min = kd_node_min(node);
    double *max = kd_node_max(node);
    double *cm = kd_node_max(node);
    double l = 0;
    double r2 = 0;

    double half[Nd];
    double full[Nd];

    if(t0->boxsize) {
        for(d = 0; d < Nd; d++) {
            half[d] = t0->boxsize[d] * 0.5;
            full[d] = t0->boxsize[d];
        }
    }

    for(d = 0; d < Nd; d++) {
        l += max[d] - min[d]; 
        double dx = pos[d] - cm[d];
        if (t0->boxsize) {
            if (dx < 0) dx = -dx;
            if (dx > half[d]) dx = full[d] - dx;
        }
        r2 += dx * dx;
    }
    double potent = 0;
    if (r2 > 0 && l*l / r2 < openangle * openangle) {
        /* no need to open the node */ 
        potent += kd_node_weight(node)[0] * gravity(sqrt(r2), userdata);
    } else {
        if(node->dim >=0) {
            potent += kd_grav(node->link[0], pos, openangle, gravity, userdata);
            potent += kd_grav(node->link[1], pos, openangle, gravity, userdata);
        } else {
            /*open the node*/ 
            double * pbase = alloca(node->size * sizeof(double) * Nd);
            double * wbase = alloca(node->size * sizeof(double) * t0->weight.dims[1]);
            ptrdiff_t i;
            double * p, * w;
            kd_collect(node, pbase, wbase);
            for(p = pbase, w = wbase, i = 0; i < node->size; i++, p+=Nd, w+= t0->weight.dims[1]) {
                r2 = 0;
                for(d = 0; d < Nd; d++) {
                    double dx = pos[d] - p[d];
                    if (t0->boxsize) {
                        if (dx < 0) dx = -dx;
                        if (dx > half[d]) dx = full[d] - dx;
                    }
                    r2 += dx * dx;
                }
                potent += w[0] * gravity(sqrt(r2), userdata);
            }

        }
    }
    return potent;
}
