#!/bin/bash

make clean

rm -f test-log

export LANG=en_US.UTF-8
export LC_ALL=en_US

make overlap >> test-log 2>&1 && \
make build_micro_hskuan >> test-log 2>&1 && \
echo Success
