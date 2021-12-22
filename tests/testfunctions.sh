export OMP_NUM_THREADS=1
export OMPI_MCA_rmaps_base_no_oversubscribe=0
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_mpi_show_mca_params=1


die () {
    exit 1
}

NUMFAILS=0

# assert command $1 returns 0.
assert_success() {
  local source=$0
  local cmd=$1
  echo "$source:" "$cmd"
  ( eval $cmd )
  if [ $? -eq 0 ]; then
     echo "[ OK   ]" "$source:" "$cmd"
     return 0
  else
     echo "[ FAIL ]" "$source:" "$cmd"
     NUMFAILS=$(($NUMFAILS + 1))
     return 1
  fi

}
# look for $2 in the file $1, die if not found.
assert_file_contains () {
  local source=$0
  local log=$1
  local pattern=$2
  if grep -q "$pattern" "$log" ; then
     echo "[ OK   ]" "$source:" Regex pattern "$pattern" matched by file "$log".
     return 0
  else
     echo "[ FAIL ]" "$source:" Regex pattern "$pattern" not matched by file "$log".
     NUMFAILS=$(($NUMFAILS + 1))
     return 1
  fi
}

report_test_status () {
  if [ $NUMFAILS -gt 0 ]; then
     echo "$0: Some tests have failed."
     die
  else
     echo "$0: All tests have passed."
  fi
}
