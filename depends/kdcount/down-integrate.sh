#
# Run this script to download MP-Sort source code.
#

tmp=$(mktemp -d)
if ! git diff --cached --exit-code; then
  echo "There are changes in the index. Commit or reset first."
  exit 1
fi

trap "rm -rf $tmp" EXIT
git clone https://github.com/rainwoodman/kdcount $tmp
sha=$(cd $tmp; git rev-parse --verify HEAD)
rsync -avz -f '- /*/' -f '- *.py' -f '- *.pyx' -f '- pykdcount.c' $tmp/kdcount/* .
git add *.[ch] Makefile
git commit -m "Update kdcount to $sha"
