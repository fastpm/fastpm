fail () {
    exit 1
}

# look for $2 in the file $1, die if not found.
expect_contains () {
  local log=$1
  local pattern=$2
  if ! grep "$pattern" "$log"; then
     echo Regex pattern "$pattern" not matched by file "$log".
     fail
  fi
  return 0
}
