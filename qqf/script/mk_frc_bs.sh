#!/bin/bash

DIR=`pwd`

kmo=12
navg=3

khpr=1
hamp=1
xdil=40.
ydil=12.
xcnt=210.
ycnt=0.

kvpr=2
vamp=8.
vdil=20.
vcnt=0.45

ovor=f
odiv=f
otmp=t
ops=f
osph=t

cd ../../solver/util/

[[ -f SETPAR ]] && rm -rf SETPAR

cat > ./SETPAR  << SETPAR 
 &nmfgt cfs='/home/qqf/Model/LBM/qqf/Gtools/psi',
        cfc='/home/qqf/Model/LBM/qqf/Gtools/chi',
        cfu='/home/qqf/Model/LBM/qqf/Gtools/u',
        cfv='/home/qqf/Model/LBM/qqf/Gtools/v',
        cfw='/home/qqf/Model/LBM/qqf/Gtools/w',
        cft='/home/qqf/Model/LBM/qqf/Gtools/t',
        cfz='/home/qqf/Model/LBM/qqf/Gtools/z',
        cfp='/home/qqf/Model/LBM/qqf/Gtools/p',
        cfq='/home/qqf/Model/LBM/qqf/Gtools/q',
        cfx='/home/qqf/Model/LBM/qqf/Gtools/dt',
        cfy='/home/qqf/Model/LBM/qqf/Gtools/dq',
        cfo='/home/qqf/Model/LBM/qqf/out/linear.t42l20.qqf.grd',
        fact=1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
        opl=t,
 &end

 &nmbs  cbs0='/home/qqf/Model/LBM/qqf/bs/qqf.t42l20',
        cbs='/home/qqf/Model/LBM/qqf/bs/qqf.t42l20.grd'
 &end

 &nmncp cncep='/home/qqf/Model/LBM/qqf/mybs/ncep.clim.y79-14.t42.grd',
        cncep2='/home/qqf/Model/LBM/qqf/mybs/ncep.clim.y79-14.ps.t42.grd',
        calt='/home/qqf/Model/LBM/bs/gt3/grz.t42',
        kmo=$kmo, navg=$navg, ozm=f, osw=f, ousez=t
 &end

 &nmecm cecm='/home/qqf/Model/LBM/bs/ecmwf/ERA40.clim.t42.grd',
        calt='/home/qqf/Model/LBM/bs/gt3/grz.t42',
        kmo=6, navg=3, ozm=f, osw=f
 &end

 &nmfin cfm='/home/qqf/Model/LBM/qqf/frc/frc.t42l20.CNP.mat',
        cfg='/home/qqf/Model/LBM/qqf/frc/frc.t42l20.CNP.grd'
        fact=1.0,1.0,1.0,1.0,1.0
 &end

 &nmvar ovor=$ovor, odiv=$odiv, otmp=$otmp, ops=$ops, osph=$osph
 &end

 &nmhpr khpr=$khpr,
        hamp=$hamp,
        xdil=$xdil,
        ydil=$ydil,
        xcnt=$xcnt,
        ycnt=$ycnt
 &end

 &nmvpr kvpr=$kvpr,
        vamp=$vamp,
        vdil=$vdil,
        vcnt=$vcnt
 &end

 &nmall owall=t
 &end

 &nmcls oclassic=t
 &end


 &nmred cdr='/home/qqf/Model/LBM/matrix.moi',
        cfo='/home/qqf/Model/LBM/matrix.moi/mat/MAT.t21l11m6.ncepannzm.moi.tmp.dat'
 &end
SETPAR

rm -rf fort.*

./ncepsbs #>& /dev/null

rm -rf fort.*

./mkfrcng #>& /dev/null

rm -rf fort.*

cd $DIR
