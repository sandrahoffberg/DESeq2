#!/bin/bash

if [ -z ${1} ]; then
    alpha="0.5"
else
    alpha=${1}
fi

if [ -z ${2} ]; then
    colNonSig="darkgray"
else
    colNonSig=${2}
fi

if [ -z ${3} ]; then
    colSig="blue"
else
    colSig=${3}
fi