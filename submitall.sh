#!/bin/bash
nohup make -j12 > makelog 2>&1 < /dev/null &
