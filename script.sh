#!/bin/bash
h5topng -t 0:166 -R -Zc dkbluered -a yarg -A nauka-eps-000000.00.h5 nauka-ezSim1.h5

convert *.png ezSim1.gif

rm *.png

h5topng -t 0:166 -R -Zc dkbluered -a yarg -A nauka-eps-000000.00.h5 nauka-ezSim2.h5

convert *.png ezSim2.gif

rm *.png

rm *.h5