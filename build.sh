#!/bin/bash
set -xe

cc -o nanomux nanomux.c tpool.c xmem.c cpthread.c -lpthread -lz -O3
cc -o nanotrim nanotrim.c tpool.c xmem.c cpthread.c -lpthread -lz -O3

# build windows
gcc -o nanomux nanomux.c tpool.c xmem.c -lpthread -lz -L. -O3
gcc -o nanotrim nanotrim.c tpool.c xmem.c -lpthread -lz -L. -O3


