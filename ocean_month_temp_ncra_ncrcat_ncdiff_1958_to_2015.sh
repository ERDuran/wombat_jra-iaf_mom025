#!/bin/bash

I_path=/g/data1a/v45/mtc599/mom5/dec16b/OUTPUTr0/

O_path=/g/data/e14/erd561/wombat_jra_mom025/

file_name1=ocean_month_
file_name2=_01
file_name2_bis=_07


years_1958_to_2014=($(seq -w 1958 1 2014))

data_path=()
idx=0
for year in ${years_1958_to_2014[@]}; do
	year_path=()
	if [ $year == 1984 ] || [ $year == 1985 ] || [ $year == 1986 ]; then
		year_path=($I_path$file_name1$year$file_name2.nc $I_path$file_name1$year$file_name2_bis.nc)
	else
		year_path=$I_path$file_name1$year$file_name2.nc
	fi

	#ncra -O -d yt_ocean,-50.0,-10.0 -d xt_ocean,-260.0,-190.0 -v temp ${year_path[@]} ${O_path}${file_name1}temp_ncra_$year.nc

	data_path[$idx]=${O_path}${file_name1}temp_ncra_$year.nc
	idx=$(expr $idx + 1)
	echo $year OK !
done

#ncrcat -O ${data_path[@]} ${O_path}${file_name1}temp_ncrcat_1958_to_2014.nc

ncwa -O -a Time ${O_path}${file_name1}temp_ncrcat_1958_to_2014.nc ${O_path}${file_name1}temp_ncwa_1958_to_2014.nc

ncdiff -O ${O_path}${file_name1}temp_ncrcat_1958_to_2014.nc ${O_path}${file_name1}temp_ncwa_1958_to_2014.nc ${O_path}${file_name1}temp_ncrcat_1958_to_2014_ncdiff_ncwa.nc
