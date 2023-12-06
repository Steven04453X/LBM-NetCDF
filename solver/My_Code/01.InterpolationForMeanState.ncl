; =======================================================
; Author: Ziming Chen 
; This code is for:
;   - interpolate the input data onto specific resolution 
;   1. indicate the input and outpur directory, and the reference data 
;   2. interpolation 
;   3. output 
; =======================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
begin
;
; =======================================================
;   1. indicate the input and outpur directory, and the reference data 
s_DirRead           = "./"
s_FileRead          = "mean_state.nc"
FileInput           = addfile(s_DirRead + s_FileRead, "r")
;
s_DirOut            = "$LNHOME/bs/test/"
s_FileOut           = s_DirOut + s_FileRead
system("rm -rf " + s_FileOut)
FileOut             = addfile(s_FileOut, "c")
GAttr               = True
copy_VarAtts(FileInput, GAttr)
GAttr@historiy      = systemfunc("date") + "  " + systemfunc("pwd") + \
                    "/01.InterpolationForMeanState.ncl"
fileattdef(FileOut, GAttr)
;
s_VarNameAll        = (/ "z", "t", "u", "v", "q", "omg", "rh" /)
;
;
s_FileRef           = "./CMIP6.historical.ltm.y1995-2014.t42.Ziming_2Models.nc"
FileRef             = addfile(s_FileRef, "r")
;
; =======================================================
;   2. interpolation 
;   3. output 
opt                 = True
opt@DstGridLon      = FileRef->lon ;r_lonMPI
opt@DstGridLat      = FileRef->lat ;r_latMPI
opt@DstReginal      = True
opt@InterpMethod    = "conserve"
opt@ForceOverwrite  = True

do iVar = 0, dimsizes(s_VarNameAll)-1
    s_VarName2      = s_VarNameAll(iVar)
    r_Var          := FileInput->$s_VarName2$
    ;
    opt@SrcGridLon := r_Var&lon
    opt@SrcGridLat := r_Var&lat 
    ;
    r_Var          := ESMF_regrid(r_Var, opt)
    system("rm -rf ./destination_grid_file.nc ./PET0.RegridWeightGen.Log ./source_grid_file.nc ./weights_file.nc")
    ;
    FileOut->$s_VarName2$ = r_Var
end do
;
print(s_FileOut + "")
;
end