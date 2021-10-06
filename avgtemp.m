
function [f] = avgtemp(albedo,emissivity,absorptivity,dissipation)
A0 = 0.02; % [m^2] Cross Sectional Area
sigma = 5.6703e-8; % [W/m^2*K^4] Stefan Constant
solar_flux = 1353; % [W/m^2] Solar flux constant approximated
r=(A0/pi)^.5;
R=6378140; % Earth radius [m]
h=300000; % Satellite altitude [m]
correction = 0.717;
Qsun = solar_flux*A0*absorptivity; % solar power directly from Sun [W]
H=h+R;
Et=((1-albedo)/4)*solar_flux;
Qir= 2*pi*r^2*Et*(1-((2*R*h + h^2)^(.5)/H)); % radiant infrared power emitted from the Earth [W]
Qearth=1.75*absorptivity; %Intperpolated approximately from Earth Reflected Solar Radiation Input to Spherical Satellites @ theta_s = 90 degrees [W]
f1= ((Qsun+Qearth+dissipation)/(emissivity)+Qir)/(sigma*A0);
f = correction*(f1)^(1/4); % [Degrees Kelvin]
end


% thm = acos(1/H)
% ths = pi/2;
%th= pi:0.01745329:2*pi;
%for i = 1:numel(th)
% syms th
% A = (2/3)*cos(ths)*(((1+(1/H^2))*(-2*H^2+1)+(2*H-1/H)*sin(ths)+(sin(ths))^2)/(H^2+1-2*H*sin(ths)^.5 + (1/H^2)+2*H))
% B = (2/pi)*vpaintegral(((H*cos(th)-1)*cos(ths)*sin(th)*cos(th))*acos(-(cos(th)*cos(ths))/(sin(th)*sin(ths)))/(H^2+1-2*H*cos(th))^1.5,th,[0 thm])
% C = (2/pi)*vpaintegral(((H^2*cos(th)-1)*cos(ths)*sin(th)*cos(th))/(H^2+1-2*H*cos(th))^1.5,th,[0 thm]);
% Qearth= solar_flux*albedo*(A+B+C);
%end

% fun4 = @(th)
% D = (2/pi)*integral(fun4,ths-pi/2,thm);
% fun5 = @(th)
% E = (2/pi)*integral(fun5,ths-pi/2,thm);
% F = (2/3)*[(2*H + 1/H^2) - (2+ 1/H^2)*(H^2-1)^.5]*cos(ths);
% if ths >= 0 && ths <=(pi/2)-thm
% Qearth = solar_flux*albedo*F
% elseif ths > ((pi/2)-thm) && ths<=(pi/2)
   
% elseif ths > pi/2 && ths<=((pi/2)+thm)
% Qearth = solar_flux*albedo*(D+E)
% end






