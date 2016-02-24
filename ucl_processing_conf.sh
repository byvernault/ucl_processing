# UCL Processing with dax configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent

UCL_PROCESSING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
UCL_PROCESSORS_DIR=$UCL_PROCESSING_DIR/ucl_processors
UCL_MODULES_DIR=$UCL_PROCESSING_DIR/ucl_modules
SPIDERS_DIR=$UCL_PROCESSING_DIR/ucl_spiders

# Edit the PATH and PYTHONPATH
export PATH=$SPIDERS_DIR:$PATH
export PYTHONPATH=$UCL_PROCESSORS_DIR:$UCL_MODULES_DIR:$PYTHONPATH

