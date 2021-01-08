#!/bin/bash
# run from parent of project directory named "dispersion-toolkit"

version=$1
dir=$2

tar --exclude=".*" -czvf dptk-$version.tgz $dir
zip -x "*.git*" -x "*.vscode*" -x ".gitignore" -r dptk-$version.zip $dir
