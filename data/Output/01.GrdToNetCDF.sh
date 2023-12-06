#!/bin/bash
#
FileInput=./TestRun.grd
#
FileForcing=frc.t42l20.Ziming.Heating_CP20N150W.nc
FileOutput=./Output_${FileForcing}
#
cdo -f nc import_binary TestRun.ctl ${FileOutput}
