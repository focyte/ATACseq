#!/bin/bash

find . -type f | egrep "[.]fq.gz" >> sample-list.txt

echo "Sample list generated successfully: sample-list.txt"
