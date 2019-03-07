load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

if (.not. isdefined("file_prefix")) then
  file_prefix = systemfunc("a=$(ls *.h0.*.nc | head -n1) && echo ${a/\.h0.*/}")
  system("echo '[Notice]: Use file_prefix " + file_prefix + ".'")
end if

file_names = systemfunc("ls " + file_prefix + "*h0*.nc")

f = addfile(file_names(0), "r")

wks = gsn_open_wks("pdf", file_prefix)

res_h = True
res_h@gsnDraw = False
res_h@gsnFrame = False
res_h@cnLinesOn = False
res_h@cnFillOn = True
res_h@cnLevelSelectionMode = "ManualLevels"
res_h@cnMinLevelValF = 5050
res_h@cnMaxLevelValF = 5950
res_h@cnLevelSpacingF = 50
; res_h@gsnPolar = "NH"
res_h@mpOutlineOn = False
res_h@gsnLeftString = "Height (m)"
res_h@sfXArray = f->lonCell
res_h@sfYArray = f->latCell

res_u = True
res_u@gsnDraw = False
res_u@gsnFrame = False
res_u@cnLinesOn = False
res_u@cnFillOn = True
res_u@cnLevelSelectionMode = "ManualLevels"
res_u@cnMinLevelValF = -20
res_u@cnMaxLevelValF =  20
res_u@cnLevelSpacingF = 2
res_u@gsnPolar = "NH"
res_u@gsnLeftString = "Edge normal wind speed (m/s)"
res_u@sfXArray = f->lonEdge
res_u@sfYArray = f->latEdge

res_panel = True
; res_panel@gsnMaximize = True
res_panel@gsnPanelRowSpec = True

do t = 0, dimsizes(file_names) - 1
  system("echo '=> " + file_names(t) + "'")
  f = addfile(file_names(t), "r")
  res_h@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  plot0 = gsn_csm_contour_map(wks, f->h(0,:), res_h)

  ; res_u@gsnRightString = cd_string(f->time, "Day %d Hour %H Min %M")
  ; plot1 = gsn_csm_contour_map(wks, f->u(0,:), res_u)

  gsn_panel(wks, (/plot0/), (/1,1/), res_panel)
end do
