tmp=$(mktemp -d)
trap "rm -rf $tmp" EXIT

git clone https://github.com/rainwoodman/MP-Sort $tmp
rsync -avz $tmp/* .
