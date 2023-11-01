#!/bin/bash

# https://medium.com/@janloo/github-actions-detected-dubious-ownership-in-repository-at-github-workspace-how-to-fix-b9cc127d4c04
sh -c "git config --global --add safe.directory $PWD"

echo pwd
pwd

make $1

echo Complete
