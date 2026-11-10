# SPDX-License-Identifier: LGPL-3.0-or-later
# solid_xl: the EIC software image (eic_xl) + the SoLID geometry (solid_dd4hep), optional Gaudi/k4FWCore.
#
#   docker build -f containers/Dockerfile -t solid_xl:dev .                       # from the repository root
#   docker build -f containers/Dockerfile --build-arg EIC_TAG=26.11.0-stable -t solid_xl:dev .
#   docker build -f containers/Dockerfile --build-arg WITH_GAUDI=1 -t solid_xl:dev-gaudi .
#   apptainer build solid_xl.sif docker-daemon://solid_xl:dev                        # for ifarm / laptops
#
# Published by .github/workflows/image.yml as ghcr.io/<repository owner>/solid_xl:<tag>; get it with
#   containers/install.sh --solid-xl
#
# solid_recon (EICrecon plugins) will be added here once it is re-based out of tree (docs/reconstruction.md).
ARG EIC_TAG=26.10.0-stable
FROM eicweb/eic_xl:${EIC_TAG}
ARG WITH_GAUDI=0
SHELL ["/bin/bash", "-c"]

# (optional) Gaudi + k4FWCore in a separate spack environment/view, reusing the image's packages
# (verified to concretize to gaudi@39.4 + k4fwcore@1.5 with 12 new packages on eic_xl 26.10).
COPY containers/spack-gaudi.yaml /opt/spack-environment/solid-gaudi/spack.yaml
RUN if [ "${WITH_GAUDI}" = "1" ]; then \
      source /opt/spack/share/spack/setup-env.sh && \
      spack -e /opt/spack-environment/solid-gaudi install --fail-fast && \
      spack -e /opt/spack-environment/solid-gaudi view --dependencies no symlink /opt/k4 gaudi k4fwcore && \
      spack clean -a ; \
    fi

# Tools eic_xl does not ship:
#   ffmpeg                       MP4 output of the event-display animations
#   gmsh (+ two X11 libraries)   mesher of the field solver, magnet_field/solve_rz.py
RUN apt-get update && \
    apt-get install -y --no-install-recommends ffmpeg libxcursor1 libxinerama1 && \
    rm -rf /var/lib/apt/lists/* && \
    python3 -m pip install --no-cache-dir --no-deps "gmsh>=4.8" && \
    ffmpeg -version | head -1 && \
    python3 -c "import gmsh; gmsh.initialize(); print('gmsh', gmsh.__version__); gmsh.finalize()"

# SoLID geometry, built with the image's own toolchain
COPY . /opt/src/solid_dd4hep
RUN cmake -S /opt/src/solid_dd4hep -B /tmp/build -DCMAKE_INSTALL_PREFIX=/opt/solid \
          -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF && \
    cmake --build /tmp/build -j"$(nproc)" && cmake --install /tmp/build && rm -rf /tmp/build

COPY containers/profile.d/z20_solid.sh /etc/profile.d/z20_solid.sh
ARG IMAGE_SOURCE=https://github.com/JeffersonLab/solid_dd4hep
ARG EIC_TAG
LABEL org.opencontainers.image.source="${IMAGE_SOURCE}" \
      org.opencontainers.image.title="SoLID runtime image (eic_xl + solid_dd4hep)" \
      org.opencontainers.image.base.name="docker.io/eicweb/eic_xl:${EIC_TAG}"
