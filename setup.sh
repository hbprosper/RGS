export RGS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=$RGS_PATH/python:$PYTHONPATH
export LD_LIBRARY_PATH=$RGS_PATH/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$RGS_PATH/lib:$DYLD_LIBRARY_PATH


