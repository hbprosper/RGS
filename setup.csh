setenv RGS_PATH `pwd`
setenv PATH ${RGS_PATH}/bin:${PATH}
setenv PYTHONPATH ${RGS_PATH}/python:${PYTHONPATH}
setenv LD_LIBRARY_PATH ${RGS_PATH}/lib:${LD_LIBRARY_PATH}
echo $RGS_PATH
