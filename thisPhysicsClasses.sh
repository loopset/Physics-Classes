#!/bin/sh

# Self locate script when sourced
if [ -n "$BASH_VERSION" ]; then
    SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
elif [ -n "$ZSH_VERSION" ]; then
    SCRIPT_DIR=$( cd -- "$( dirname -- "${0}" )" &> /dev/null && pwd )
else
    echo "thisPhysicsClasses.sh: unsupported shell to source"
    exit 1
fi

# Add to system libraries
lib="${SCRIPT_DIR}/install/lib"
if [[ "$(uname)" == "Darwin" ]]; then # MacOS
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${lib}"
elif [[ "$(uname)" == "Linux" ]]; then # GNU/Linux
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${lib}"
else
    echo "thisPhysicsClasses.sh: unsupported system to source"
    exit 1
fi
# And lastly to ROOT includes
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:"${SCRIPT_DIR}/install/include"
