#!/bin/bash
set -xe

cc -o nanomux thpool.c nanomux.c -lz -lpthread -O3 
cc -o nanotrim nanotrim.c thpool.c -lz -lpthread -lm -O3

