#!/bin/bash -e
shopt -s nullglob
unalias -a

build_library() {
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
   comp_flags=("${base_flags[@]}")
   if $debug_build; then
      comp_flags+=("${build_debug_flags[@]}")
   else
      comp_flags+=("${build_flags[@]}")
   fi
   while read -r filename; do
      srcfile=$srcdir/$filename
      objfile=$blddir/${filename%.*}.o
      if ! test -e "$objfile" || test "$srcfile" -nt "$objfile"
      then
         echo Compiling $filename...
         "$FC" "${comp_flags[@]}" -J "$blddir" -c "$srcfile" -o "$objfile"
      fi
      object_files+=("$objfile")
   done < <(grep -Ehv '^$|^#' "$srcdir/program_dependencies.txt")
}

build_programs() {
   while read -r filename progname; do
      srcfile=$srcdir/$filename
      execfile=$blddir/$progname
      echo Building program $progname...
      "$FC" "$srcfile" "${object_files[@]}" -o "$execfile" "${comp_flags[@]}" "${link_flags[@]}" -J "$blddir"
   done < <(grep -Ehv '^$|^#' "$srcdir/program_list.txt")
}

blddir=$PWD/build
srcdir=$PWD/fortran

if test ! -e "$blddir"; then
   mkdir "$blddir"
elif test ! -d "$blddir"; then
   echo Error: $blddir exists but is not a directory
   exit 1
fi

## Compile with gfortran
FC=gfortran
base_flags=(-std=f2008)
build_flags=(-O3)
build_debug_flags=(-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -Wunused)

## Compile with ifort
#FC=ifort
#base_flags=(-stand f08)
#build_flags=(-O3)
#build_debug_flags=(-O0 -g -traceback -fpe0 -check bounds -warn all)

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
