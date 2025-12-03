FROM ubuntu:20.04 AS build

RUN apt-get update && apt-get install -y g++ libz-dev

# Copy dependencies
COPY ccp4srs /usr/src/pisa/ccp4srs
COPY mmdb2 /usr/src/pisa/mmdb2
COPY pisa /usr/src/pisa/pisa
COPY pisalib /usr/src/pisa/pisalib
COPY ssm /usr/src/pisa/ssm
COPY compile.sh /usr/src/pisa/compile.sh

WORKDIR /usr/src/pisa

ENV SRCDIR=/usr/src/pisa

# Compile
RUN chmod +x /usr/src/pisa/compile.sh
RUN /usr/src/pisa/compile.sh

FROM ubuntu:20.04

# Copy binary from build
COPY --from=build /usr/src/pisa/build /usr/bin/pisa

# Copy config files
COPY setup /usr/share/pisa/setup
ENV PISA_CONFIG_DIR=/usr/share/pisa/setup

# Update constants in main config file
ENV DATA_DIR=/data
RUN mkdir -p ${DATA_DIR}

RUN sed -i "s|path_dataroot|${DATA_DIR}|g" ${PISA_CONFIG_DIR}/pisa_cfg_tmp && \
    sed -i "s|path_to_setup/srs|${PISA_CONFIG_DIR}/srs|g" ${PISA_CONFIG_DIR}/pisa_cfg_tmp && \
    sed -i "s|path_to_setup/molref|${PISA_CONFIG_DIR}/molref|g" ${PISA_CONFIG_DIR}/pisa_cfg_tmp && \
    sed -i "s|path_to_setup/pisastore|${PISA_CONFIG_DIR}/pisastore|g" ${PISA_CONFIG_DIR}/pisa_cfg_tmp && \
    sed -i "/pisa_/d" ${PISA_CONFIG_DIR}/pisa_cfg_tmp

ENV PATH=$PATH:/usr/bin/pisa
