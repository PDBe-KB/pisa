FROM ubuntu:20.04 AS build

RUN apt-get update && apt-get install -y g++ libz-dev

# Copy dependencies
COPY ccp4srs /usr/src/pisa-lite/ccp4srs
COPY mmdb2 /usr/src/pisa-lite/mmdb2
COPY pisa /usr/src/pisa-lite/pisa
COPY pisalib /usr/src/pisa-lite/pisalib
COPY ssm /usr/src/pisa-lite/ssm
COPY compile.sh /usr/src/pisa-lite/compile.sh

WORKDIR /usr/src/pisa-lite

ENV SRCDIR=/usr/src/pisa-lite

# Compile
RUN chmod +x /usr/src/pisa-lite/compile.sh
RUN /usr/src/pisa-lite/compile.sh

FROM ubuntu:20.04

COPY --from=build /usr/src/pisa-lite/build /usr/share/pisa-lite

ENV PATH=$PATH:/usr/share/pisa-lite
