#!/bin/bash

I_path=/g/data1a/v45/mtc599/mom5/dec16b/OUTPUTr0/

O_path=/g/data/e14/erd561/wombat_jra_mom025/

file_name1=ocean_surface_
file_name2=_01

dec0614=($(seq -w 2006 1 2014))
dec9705=($(seq -w 1997 1 2005))
dec8896=($(seq -w 1988 1 1996))
dec7685=($(seq -w 1976 1 1985))
dec6675=($(seq -w 1966 1 1975))

for f in {0..9}; do
	dec0614_paths[${f}]=${I_path}${file_name1}${dec0614[${f}]}${file_name2}.nc
	dec9705_paths[${f}]=${I_path}${file_name1}${dec9705[${f}]}${file_name2}.nc
	dec8896_paths[${f}]=${I_path}${file_name1}${dec8896[${f}]}${file_name2}.nc
done

ncra -v temp ${dec0614_paths} ${O_path}${file_name1}temp_ncra2006to2014.nc
ncra -v temp ${dec9705_paths} ${O_path}${file_name1}temp_ncra1997to2005.nc
ncra -v temp ${dec8896_paths} ${O_path}${file_name1}temp_ncra1988to1996.nc

ncdiff ${O_path}${file_name1}temp_ncra2006to2014.nc ${O_path}${file_name1}temp_ncra1997to2005.nc ${O_path}${file_name1}temp_diff_ncra2006to2014_ncra1997to2005.nc
ncdiff ${O_path}${file_name1}temp_ncra1997to2005.nc ${O_path}${file_name1}temp_ncra1988to1996.nc ${O_path}${file_name1}temp_diff_ncra1997to2005_ncra1988to1996.nc

