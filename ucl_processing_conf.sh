# UCL Processing with dax configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent


UCL_PROCESSING_DIR=/home/byvernau/Xnat-management/ucl_processing
SPIDERS_DIR=$UCL_PROCESSING_DIR/spiders

# Edit the PATH and PYTHONPATH
export PATH=$SPIDERS_DIR:$PATH

export PYTHONPATH=$UCL_PROCESSING_DIR:$PYTHONPATH

