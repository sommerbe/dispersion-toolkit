#!/bin/bash
# run from <project root dir>

# audience:
# intended for source developers, NOT FOR BUILDERS / USERS
# i.e. developers are those among us who modify manpage.md files
#

# dependencies:
# pandoc, man, python

script_configure_file="./script/configure-file.py"

function mk_manpage_pre() {
  # args
  pdir=$1
  section=$2

  # source files
  src_md="$pdir/manpage.md"
  src_tpl="$pdir/manpage.hpp.in"

  # intermediate files
  tmp_pre="$pdir/manpage.$section.pre"

  # final output
  dst_man="$pdir/manpage.$section"
  dst_hpp="$pdir/manpage.hpp"

  # convert manpage mardown to man roff
  pandoc -s -t man "$src_md" -o "$dst_man"

  # generate pre-format to be embedded into binary --help output
  man -l "$dst_man" > "$tmp_pre"

  # generate include header
  python "$script_configure_file" "$tmp_pre" "$src_tpl" manpage > "$dst_hpp"

  # cleanup intermediate files
  rm "$tmp_pre"

  # log
  echo "mkman $pdir"
}

# check: correct pwd
if [ ! -f "$script_configure_file" ]; then
  echo "fatal: consider running this script from project's root directory."
  exit 1
fi

mk_manpage_pre "src/set/fibonacci" 1
mk_manpage_pre "src/set/kritzinger" 1
mk_manpage_pre "src/set/cswap" 1
mk_manpage_pre "src/measure/dispnaamad" 1
mk_manpage_pre "src/measure/dispdumitrescu2017" 1
mk_manpage_pre "src/measure/dispcombinatorial" 1
mk_manpage_pre "src/measure/dispgs" 1
mk_manpage_pre "src/measure/pdisppermute" 1
mk_manpage_pre "src/opt/mindispgs" 1
mk_manpage_pre "src/stat/confidence" 1
mk_manpage_pre "src/vis/psspy" 1
mk_manpage_pre "src/adapter/utk" 1
mk_manpage_pre "src/adapter/readmatrix" 1
mk_manpage_pre "src/adapter/writematrix" 1
