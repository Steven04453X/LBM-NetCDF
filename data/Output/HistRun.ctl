* products of historical run
* steady response to zero forcing 
* gamma=6000, solution for EPS=1.D-2
DSET ^HistRun.grd
OPTIONS SEQUENTIAL YREV
TITLE time-integration
UNDEF -999.
*UNDEF -999999.
OPTIONS big_endian
XDEF 128 LINEAR 0. 2.81250
YDEF 64  LEVELS 
-87.864 -85.097 -82.313 -79.526 -76.737 -73.948 -71.158 -68.368 -65.578 
-62.787 -59.997 -57.207 -54.416 -51.626 -48.835 -46.045 -43.254 -40.464 
-37.673 -34.883 -32.092 -29.301 -26.511 -23.720 -20.930 -18.139 -15.348 
-12.558  -9.767  -6.976  -4.186  -1.395   1.395   4.186  6.976   9.767  
12.558  15.348  18.139  20.930  23.720  26.511  29.301 32.092  34.883  
37.673  40.464  43.254  46.045  48.835  51.626  54.416  57.207  59.997  
62.787  65.578  68.368  71.158  73.948  76.737  79.526  82.313  85.097  
87.864 
ZDEF 20  LEVELS 1000 950 900 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5
TDEF 29 LINEAR 1jan0000 1dy
VARS 8
psi    20 99 stream function     [m**2/s]
chi    20 99 velocity potential  [m**2/s]
u      20 99 zonal wind          [m/s]
v      20 99 meridional wind     [m/s]
w      20 99 p-vertical velocity [hPa/s]
t      20 99 temperature         [K]
z      20 99 geopotential height [m]
p       1 99 surface pressure    [hPa]
ENDVARS
