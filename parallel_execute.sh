#!/bin/bash
cd $1
parallel ../run.sh ::: harmonic-*
cd ..
