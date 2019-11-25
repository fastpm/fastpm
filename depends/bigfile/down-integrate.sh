#
# Run this script to download MP-Sort source code.
#

tmp=$(mktemp -d)
if ! git diff --cached --exit-code; then
  echo "There are changes in the index. Commit or reset first."
  exit 1
fi

trap "rm -rf $tmp" EXIT
git clone https://github.com/rainwoodman/bigfile $tmp
sha=$(cd $tmp; git rev-parse --verify HEAD)
rsync -avz -f '- /*/' $tmp/src/* .
git add *.[ch] Makefile *.txt
git commit -m "Update bigfile to $sha"
