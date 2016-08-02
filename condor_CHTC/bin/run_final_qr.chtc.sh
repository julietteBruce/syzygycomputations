#!/bin/sh
# script for execution of deployed applications
#


## Set up local matlab environment;
wget http://proxy.chtc.wisc.edu/SQUID/r2014b.tar.gz;
tar zxf r2014b.tar.gz;

## Untar the data matrix;
dataFile=$2
outFile=$3

/bin/gunzip $dataFile;
dataFile=${dataFile/\.gz/};




##### This environment set up is the output of mcc.

# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#


exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=

  while [ $# -gt 0 ]; do
      token=$1
      args="${args} \"${token}\"" 
      shift
  done
  echo "Executing \"${exe_dir}/final_qr\"" $dataFile $outFile
  eval "\"${exe_dir}/final_qr\"" $dataFile $outFile
fi
exit

