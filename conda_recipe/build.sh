#!/bin/bash

## change to source dir
cd ${SRC_DIR}

## compile
make

## install
mkdir -p $PREFIX/bin
cp -r ${SRC_DIR}/* ${PREFIX}/bin/
