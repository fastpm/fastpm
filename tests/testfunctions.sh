fail () {
    exit 1
}
export OMP_NUM_THREADS=1

FASTPM=`dirname $0`/../src/fastpm
set -x
