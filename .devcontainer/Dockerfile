# Attempt to replicate github action runner environment
FROM ghcr.io/catthehacker/ubuntu:runner-20.04 AS binary_builder


#Dependencies for binary build
RUN sudo -- sh -c 'apt update && export DEBIAN_FRONTEND=noninteractive \
     && apt -y install --no-install-recommends libwxgtk3.0-gtk3-dev libproj-dev libavcodec-dev'

# Additional dependencies for source distribution build
FROM binary_builder AS source_dist_builder

RUN sudo -- sh -c 'add-apt-repository universe && \
                  add-apt-repository ppa:inkscape.dev/stable && \
                  apt-get update && export DEBIAN_FRONTEND=noninteractive && \
                  apt-get -y install --no-install-recommends inkscape unifont docbook-utils netpbm w3m liblocale-po-perl'