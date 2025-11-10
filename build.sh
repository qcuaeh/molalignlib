#!/bin/bash -e
shopt -s nullglob
unalias -a

build_library() {
   srcdir=$PWD/modules
   if $full_build; then
      full_build=false
      if test -d "$blddir"; then
         for file in "$blddir"/*{.mod,.o,.a}; do
            rm "$file"
         done
      fi
   fi
   if test ! -d "$srcdir"; then
      echo Error: $srcdir does not exist
      exit 1
   fi
   comp_flags=("${env_base_flags[@]}")
   if $debug_build; then
      comp_flags+=("${env_build_debug_flags[@]}")
   else
      comp_flags+=("${env_build_flags[@]}")
   fi
   while read -r filename; do
      srcfile=$srcdir/$filename
      objfile=$blddir/${filename%.*}.o
      if ! test -e "$objfile" || test "$srcfile" -nt "$objfile"
      then
         echo Compiling $filename...
         "$env_FC" "${comp_flags[@]}" -J "$blddir" -c "$srcfile" -o "$objfile"
      fi
      ar rcs "$blddir/molalignlib.a" "$objfile"
   done < <(grep -v '^#' "$srcdir/module_files")
}

build_programs() {
   srcdir=$PWD/programs
   molalignlib=$blddir/molalignlib.a
   while read -r filename progname; do
      srcfile=$srcdir/$filename
      execfile=$blddir/$progname
      echo Building program $progname...
      "$env_FC" "${comp_flags[@]}" "${link_flags[@]}" -J "$blddir" "$srcfile" "$molalignlib" -o "$execfile"
   done < <(grep -v '^#' "$srcdir/program_files")
}

if test ! -e ./build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f ./build.env; then
   echo Error: build.env exists but is not a file
   exit 1
fi

blddir=$PWD/build

if test ! -e "$blddir"; then
   mkdir "$blddir"
elif test ! -d "$blddir"; then
   echo Error: $blddir exists but is not a directory
   exit 1
fi

# Read compiler flags
while IFS= read -r line; do
   lhs=${line%%=*}
   rhs=${line#*=}
   IFS=\  read -r -a "env_$lhs" <<< "$rhs"
done < <(grep -v '^#' ./build.env)

#echo "$env_FC"
#echo "${env_base_flags[@]}"
#echo "${env_build_flags[@]}"
#echo "${env_build_debug_flags[@]}"

# Get build options
full_build=true
debug_build=false
while getopts ":dq" opt; do
  case $opt in
    d)
      debug_build=true
      ;;
    q)
      full_build=false
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Build the static library
build_library

# Build the RMSD programs
build_programs
