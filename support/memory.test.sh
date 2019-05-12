#!/bin/bash
$@ &
pid=$!
echo "Command: $@" >> ./memory_usage_$pid
trap "kill $pid 2> /dev/null" EXIT
while kill -0 $pid 2> /dev/null; do
    cat /proc/$pid/status | grep "VmHWM" | awk {'print $2'} >> ./memory_usage_$pid
    sleep 1
done
trap - EXIT
echo $(cat ./memory_usage_$pid)
