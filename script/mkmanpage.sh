#!/bin/bash
# run from <project root dir>

function mk_manpage_pre() {
  # args
  pdir=$1
  section=$2

  # files
  src_md="$pdir/manpage.md"
  src_tpl="$pdir/manpage.hpp.in"
  dst_man="$pdir/manpage.$section"
  dst_pre="$pdir/manpage.$section.pre"
  dst_hpp="$pdir/manpage.hpp"

  # convert manpage mardown to man roff
  pandoc -s -t man "$src_md" -o "$dst_man"

  # generate pre-format to be embedded into binary --help output
  man -l "$dst_man" > "$dst_pre"

  # generate include header
  python script/configure-file.py "$dst_pre" "$src_tpl" manpage > "$dst_hpp"
}

mk_manpage_pre "src/set/fibonacci" 1