#!/bin/bash

rm -rf ncepsum.t42l20.nc
cdo -f nc import_binary ../ncepsum.t42l20.ctl ncepsum.t42l20.nc
