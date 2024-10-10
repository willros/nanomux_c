#!/bin/bash
set -xe

cc -o nanomux nanomux.c tpool.c xmem.c cpthread.c -lpthread -lz -O3
cc -o nanotrim nanotrim.c tpool.c xmem.c cpthread.c -lpthread -lz -O3

