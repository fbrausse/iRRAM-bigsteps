#!/bin/bash

if [ $# -ne 1 ]; then
	echo "usage: $0 <alpha>" >&2
	exit 1
fi

alpha=$1

echo -n "\
2
0 0 0 1  1.0
1 0 0 1  $alpha
1 0 2 1 -$alpha
1 0 1 0 -1.0
"
