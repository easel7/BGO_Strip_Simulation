#!/bin/bash
Tag=Carbon
particle=ion
# Loop through energy values from 10 to 100 GeV with a step of 10
for energy in $(seq 10 10 100)
do
  # Generate the corresponding macro file for the current energy value
  cat <<EOL > "${Tag}_${energy}GeV.mac"
/run/initialize
/analysis/setFileName ../Root/${Tag}_${energy}GeV
/gps/particle ${particle}
/gps/ion 6 12 6
/gps/ene/type Mono
/gps/ene/mono ${energy} GeV

/gps/pos/type Plane
/gps/pos/shape Square
/gps/pos/halfx 10 cm
/gps/pos/halfy 10 cm
/gps/pos/centre 0. 0. -214.2 mm

/gps/direction 0 0 1

/run/printProgress 1000
/run/beamOn 10000
EOL
done

for energy in $(seq 200 100 1000)
do
  # Generate the corresponding macro file for the current energy value
  cat <<EOL > "${Tag}_${energy}GeV.mac"
/run/initialize
/analysis/setFileName ../Root/${Tag}_${energy}GeV
/gps/particle ${particle}
/gps/ion 6 12 6
/gps/ene/type Mono
/gps/ene/mono ${energy} GeV

/gps/pos/type Plane
/gps/pos/shape Square
/gps/pos/halfx 10 cm
/gps/pos/halfy 10 cm
/gps/pos/centre 0. 0. -214.2 mm

/gps/direction 0 0 1

/run/printProgress 1000
/run/beamOn 10000
EOL
done

for energy in $(seq 2 1 10)
do
  # Generate the corresponding macro file for the current energy value
  cat <<EOL > "${Tag}_${energy}TeV.mac"
/run/initialize
/analysis/setFileName ../Root/${Tag}_${energy}TeV
/gps/particle ${particle}
/gps/ion 6 12 6
/gps/ene/type Mono
/gps/ene/mono ${energy} TeV

/gps/pos/type Plane
/gps/pos/shape Square
/gps/pos/halfx 10 cm
/gps/pos/halfy 10 cm
/gps/pos/centre 0. 0. -214.2 mm

/gps/direction 0 0 1

/run/printProgress 1000
/run/beamOn 10000
EOL
done

