#!/bin/bash

while read dirname
do
	rm -r  $dirname
done < dirnames.txt
rm dirnames.txt


