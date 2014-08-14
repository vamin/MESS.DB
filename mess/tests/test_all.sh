#!/bin/bash

cd "$(dirname "$0")"
python2.7 -m unittest discover
