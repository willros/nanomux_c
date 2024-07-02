#!/bin/bash
set -xe

if [[ -z $1 ]];
then 
    echo "building WITHOUT read splitting";
    cc -o nanomux nanomux.c -lz -lpthread -O3 -Wall -Wextra -pedantic
else
    echo "building WITH read splitting";
    cc -o nanomux nanomux.c -DREAD_SPLITTING -lz -lpthread -O3 -Wall -Wextra -pedantic
fi

