DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=$DIR/python:$PYTHONPATH
export LD_LIBRARY_PATH=$DIR/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DIR/lib:$DYLD_LIBRARY_PATH
export RGS_PATH=$DIR

