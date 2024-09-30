#!/bin/bash
set -xe

cc -o nanomux nanomux.c -lz -lpthread -O3 
cc -o trim trim.c thpool.c -lz -lpthread -lm -O3

