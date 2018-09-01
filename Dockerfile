FROM ubuntu:18.04

RUN apt-get update -y && \
    apt-get install -y \
        build-essential \
        cmake \
        git \
        python3 \
        time \
        zlib1g-dev

ADD Big-BWT /code/Big-BWT
ADD CMakeLists.txt /code/CMakeLists.txt
ADD internal /code/internal
ADD klib /code/klib
ADD *.cpp /code/
ADD scripts /code/scripts
ADD sdsl-lite /code/sdsl-lite

RUN cd /code/sdsl-lite && \
        cmake . && \
        make

RUN cd /code/Big-BWT && make

RUN mkdir -p /code/release /code/debug
RUN cd /code/release && \
        cmake -DCMAKE_BUILD_TYPE=Release .. && \
        make -j2 ri-align ri-buildfasta

RUN cd /code/debug && \
        cmake -DCMAKE_BUILD_TYPE=Debug .. && \
        make -j2 ri-align ri-buildfasta
