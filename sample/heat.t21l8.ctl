* sample .ctl for heating in moist LBM
* can be obtained by mksfrc.moist.csh
DSET ^heat.t21l8.grd
OPTIONS SEQUENTIAL YREV
TITLE dumy
UNDEF -999.
XDEF 64 LINEAR 0. 5.625
YDEF 32 LEVELS -85.761 -80.269 -74.745 -69.213 -63.679 -58.143 -52.607 
-47.070 -41.532 -35.995 -30.458 -24.920 -19.382 -13.844 -8.3067 -2.7689 
2.7689 8.3067 13.844 19.382 24.920 30.458 35.995 41.532 47.070 52.607 
58.143 63.679 69.213 74.745 80.269 85.761
ZDEF 8  LEVELS 0.995 0.945 0.830 0.653 0.463 0.303 0.162 0.0415
TDEF 1 LINEAR 15jan0000 1mo
VARS 4
dtc     8 99 heat source due to convection [K/s]
dqc     8 99 moisture source due to convection [kg/kg/s]
dts     8 99 heat source due to sfc. fluxes [K/s]
dqs     8 99 moisture source due to sfc. fluxes [kg/kg/s]
ENDVARS
