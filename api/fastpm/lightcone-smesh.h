typedef struct FastPMSMesh {
    FastPMLightCone * lc;

    double smoothing;
    FastPMColumnTags attributes;

    struct FastPMSMeshLayer {
        enum {
            FASTPM_SMESH_SPHERE,
            FASTPM_SMESH_PLANE,
        } type;

        union {
            struct {
                double * ra;
                double * dec;
                uint64_t * pix;
               double (* vec)[3];
            };
            struct {
                double (* xy)[2];
            };
        };

        int Nxy;
        int nside; /* healpix number of pixels; 0 for other cases. */

        double * a;
        double * z;
        int Na;

        struct FastPMSMeshLayer * next;
    } * layers;

    size_t np_upper;

    /* state about the last time range */
    struct {
        FastPMStore p[1];
        double a_f; /* the time potential was updated */
    } last;
    int started;

    /* Extensions */
    FastPMEventHandler * event_handlers;

} FastPMSMesh;

void
fastpm_smesh_init(FastPMSMesh * mesh, FastPMLightCone * lc, size_t np_upper, double smoothing);

typedef struct {
    double aemit;
    int nside;
    double distance;
} FastPMSMeshSlice;

FastPMSMeshSlice *
fastpm_smesh_get_aemit(FastPMSMesh * mesh, size_t * Na);

struct FastPMSMeshLayer *
fastpm_smesh_add_layer_plane(FastPMSMesh * mesh,
        double (*xy)[2], size_t Nxy,
        double * a, size_t Na);

struct FastPMSMeshLayer *
fastpm_smesh_add_layer_pm(FastPMSMesh * mesh,
        PM * pm, double * shift, ptrdiff_t * Nc,
        double amin, double amax);

struct FastPMSMeshLayer *
fastpm_smesh_add_layer_sphere(FastPMSMesh * mesh,
        double * ra, double * dec, uint64_t * pix, size_t Npix,
        double * a, size_t Na);

struct FastPMSMeshLayer *
fastpm_smesh_add_layer_healpix(FastPMSMesh * mesh,
       int nside,
        double * a, size_t Na, MPI_Comm comm);

void
fastpm_smesh_add_layers_healpix(FastPMSMesh * mesh,
        double surface_density, double line_density,
        double amin, double amax, int maxnside,
        MPI_Comm comm);

int
fastpm_smesh_select_active(FastPMSMesh * layer,
        double a0, double a1,
        FastPMStore * q
    );

int
fastpm_smesh_compute_potential(
        FastPMSMesh * mesh,
        PM * pm,
        FastPMPainter * painter,
        FastPMKernelType kernel,
        FastPMFloat * delta_k,
        double a_f,
        double a_n
);

void
fastpm_smesh_destroy(FastPMSMesh * mesh);
