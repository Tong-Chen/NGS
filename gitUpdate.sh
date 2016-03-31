#!/bin/bash

#set -x
set -e
set -u

if test $# -lt 2; then
	echo 1>&2 "Usage $0 file commit"
	exit 1
fi

git add $1
git commit -m "$2 FOR <$1>" 
git push -u origin master
