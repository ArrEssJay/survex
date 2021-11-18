# Attempt to replicate github action runner environment
FROM ghcr.io/catthehacker/ubuntu:runner-20.04 AS binary_builder


#Dependencies for binary build
RUN sudo -- sh -c 'apt update && export DEBIAN_FRONTEND=noninteractive \
     && apt -y install --no-install-recommends libwxgtk3.0-gtk3-dev libproj-dev libavcodec-dev'
     
# automake pkg-config gcc g++ libwxgtk3.0-gtk3-dev libavcodec-dev libavformat-dev libproj-dev             sudo cpan -T -i Locale::PO
#mesa-common-dev libglu1-mesa-dev libx11-dev libxext-dev ghostscript netpbm x11proto-core-dev liblocale-po-perl unifont docbook-utils w3m inkscape


# Additional dependencies for source distribution build
FROM binary_builder AS source_dist_builder

RUN sudo -- sh -c 'add-apt-repository universe && \
                  add-apt-repository ppa:inkscape.dev/stable && \
                  apt-get update && export DEBIAN_FRONTEND=noninteractive && \
                  apt-get -y install --no-install-recommends inkscape unifont docbook-utils netpbm'