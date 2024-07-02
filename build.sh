#!/bin/bash
set -xe

cc -o nanomux nanomux.c -lz -lpthread -O3 -Wall -Wextra -pedantic

