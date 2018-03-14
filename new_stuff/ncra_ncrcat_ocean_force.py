'''
This file creates annual means of tau in wombat_jra-iaf_mom025 and concatenates the result in one file.
Takes data from
/g/data1a/v45/mtc599/mom5/dec16b/OUTPUTr0/
and places the outputs in
/g/data/e14/erd561/wombat_jra-iaf_mom025/

Earl Duran 
created: 14-Mar-18
e.duran@unsw.edu.au
'''

import os

input_path = '/g/data1a/v45/mtc599/mom5/dec16b/OUTPUTr0/'

file_number = list(range(1958, 2015, 1))

output_path = '/g/data/e14/erd561/wombat_jra-iaf_mom025/'

pls = os.listdir(output_path)

var = ['tau_x', 'tau_y', 'sfc_hflux']
x_axis = ['xu_ocean', 'xu_ocean', 'xt_ocean']
y_axis = ['yu_ocean', 'yu_ocean', 'yt_ocean']

for n in file_number:
    input_data_path = input_path
    if n in [1984,1985,1986]:
        input_data_path += 'ocean_force_' + str(n) + '_01.nc ' + \
        input_path + 'ocean_force_' + str(n) + '_07.nc'
    else:
        input_data_path += 'ocean_force_' + str(n) + '_01.nc'
            
    for v,x,y in zip(var,x_axis,y_axis):
        output_data = v + '_' + str(n) + '.nc'
        output_data_path = output_path + output_data

        if output_data not in pls:
            os.system('ncra -d ' + y + ',-60.0,-20.0 -d ' + x + ',-260.0,-190.0 -v ' \
                + v + ' ' + input_data_path + ' ' + output_data_path)
            print(output_data_path + ' OK')


for v,x,y in zip(var,x_axis,y_axis):
    output_data = v + '_' + str(file_number[0]) + '-' + str(file_number[-1]) + '.nc'
    output_data_path = output_path + output_data
    
    if output_data not in pls:
        
        input_data_path = ''
        for n in file_number:
            input_data_path += output_path + v + '_' + str(n) + '.nc '
        
        
        os.system('ncrcat -v ' \
            + v + ' ' + input_data_path + ' ' + output_data_path)
        print(output_data_path + ' OK')
        
        
print('DONE') 
    
    
    
    
    