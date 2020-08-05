#!/bin/bash

while read dirname
do
	echo $dirname
	cat $dirname/out*
	echo ''
done < dirnames.txt



