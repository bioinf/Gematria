#!/bin/bash

FILES='src/*.c'

last_time=`stat $FILES | grep Modify`
while true; do
  current_time=`stat $FILES | grep Modify`
  if [[ "$current_time" != "$last_time" ]]; then

    # ----------------------- #
    rm -rf makegms* dist build;
    python3 setup.py build;
    python3 setup.py install;
    # ----------------------- #

    last_time=$current_time;
  fi
  sleep 1
done;
