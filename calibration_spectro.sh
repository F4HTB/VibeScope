#!/bin/bash

# Liste des fréquences à tester (exemple avec tes bandes)
freqs=(20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000)

for f in "${freqs[@]}"; do
  echo "Lecture d'une sinusoïde à ${f} Hz"
  sox -n -c 2 -t alsa plughw:CARD=loopTest1,DEV=1 synth 5 sine $f vol 0dB
done
