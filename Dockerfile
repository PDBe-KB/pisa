FROM ubuntu:20.04

RUN apt-get update && apt-get install -y g++ libz-dev

COPY . /usr/src/pisa-lite
WORKDIR /usr/src/pisa-lite

ENV SRCDIR=/usr/src/pisa-lite

RUN chmod +x /usr/src/pisa-lite/compile.sh
RUN /usr/src/pisa-lite/compile.sh

CMD ["sh", "-c", "./build/pisa docker-pisa -analyse example_data/6gve.cif pisa_cfg && ./build/pisa docker-pisa -xml interface pisa_cfg"]