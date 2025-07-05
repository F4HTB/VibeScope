#!/bin/bash

echo "Lecture d'une sinusoïde à 1000 Hz"
sox -n -t alsa -b 16 -e signed-integer -c 2 plughw:CARD=loopTest1,DEV=1 synth 60 sine 1000 vol 0dB
