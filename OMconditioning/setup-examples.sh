#!/bin/bash
# script to run the process_gpr.R script
# abspath function to find absolute path
abspath() {
  [[ $1 = /* ]] && ABSPATH="$1" || ABSPATH="$PWD/${1#./}"
}

#ScriptFolder="$(readlink -f ${BASH_SOURCE[0]})"
abspath "${BASH_SOURCE[0]}"
ScriptFolder=$ABSPATH
ScriptFolder=${ScriptFolder//" "/"\ "}
ScriptFolder=$(dirname "$ScriptFolder")"/"

cd $ScriptFolder

unzip BET.zip
unzip YFT.zip
