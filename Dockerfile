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

COPY --from=build /usr/src/pisa/build /usr/share/pisa

ENV PATH=$PATH:/usr/share/pisa
