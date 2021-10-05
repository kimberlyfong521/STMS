clear
clc

T= readtable('Satellite4_SEET_Temperature.csv');
T_ambient=273.15-20.63;
emissivity = 0.34;
flux=[];
temp = T(:,2);
satellite_temp = table2array(temp);
sigma= 5.6703e-8;
dates = ( datetime('19-Nov-2020 8:00:00'):seconds(60):datetime('22-Nov-2020 8:00:00') )';
range = isbetween(dates,'19-Nov-2020','22-Nov-2020');
num_rows=size(temp,1);
for i=1:1:num_rows

flux(i) = emissivity*sigma*((satellite_temp(i)+273.15)^4 - T_ambient^4);

end
Heat_Flux = flux.';
mean_Heat_Flux = mean(Heat_Flux);
maxheatflux = max(Heat_Flux);
minheatflux = min(Heat_Flux);
plot(dates(range),Heat_Flux(range))
text(1,22,['Average Heat Flux = ',num2str(mean_Heat_Flux),'W/m^2']);
text(1,20,['Maximum Heat Flux = ',num2str(maxheatflux),'W/m^2']);
text(1,21,['Minimum Heat Flux = ',num2str(minheatflux),'W/m^2']);

