# SoLID defaults inside the solid_xl image (runs before the image's z21_epic_main.sh, which only
# sets up ePIC when DETECTOR_PATH/DETECTOR_CONFIG are still unset).
export SOLID_PREFIX=/opt/solid
if [ -z "${DETECTOR_PATH}" ] && [ -f /opt/solid/bin/thissolid.sh ]; then
  . /opt/solid/bin/thissolid.sh "${SOLID_DEFAULT_CONFIG:-solid}"
fi
if [ -d /opt/k4 ]; then
  export CMAKE_PREFIX_PATH="/opt/k4${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"
  export LD_LIBRARY_PATH="/opt/k4/lib:${LD_LIBRARY_PATH}" PYTHONPATH="/opt/k4/python:${PYTHONPATH}" PATH="/opt/k4/bin:${PATH}"
fi
