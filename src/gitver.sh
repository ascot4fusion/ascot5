#!/bin/bash
GIT_BRANCH="$(git branch | grep \* | cut -d ' ' -f2)"
GIT_VERSION="$(git describe --abbrev=7 --dirty --always --tags)"

if [ -z "$GIT_VERSION" ]
then
    GIT_BRANCH="main"
    GIT_VERSION=$( <VERSION)
fi

echo "#ifndef GITVER_H" > gitver.h
echo "#define GITVER_H" >> gitver.h
echo "" >> gitver.h

echo "#define GIT_BRANCH" \"$GIT_BRANCH\" >> gitver.h
echo "#define GIT_VERSION" \"$GIT_VERSION\" >> gitver.h

echo "" >> gitver.h
echo "#endif" >> gitver.h
