#!/bin/sh

RUN_ARG=
if [ $1 == "ssh" ] ; then
    RUN_ARG="-it"
fi

docker run ${RUN_ARG} -v `pwd`:/wd -v $HOME/.aws:/root/.aws benlangmead/vagrant-docker $*
