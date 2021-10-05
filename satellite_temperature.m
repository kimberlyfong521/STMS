clear
clc

% VARIABLES
h=500000;                           % Satellite altitude [m]
absorptivity = .85;
emissivity = .85;
dissipation = 3;
albedo = .31;
A_l = .02;                          % Area of 2-Unit Side
A_s = .01;                          % Area of 1-Unit Side
A_tot = A_l * 4 + A_s * 2;          % Total Surface Area (m^2)
Cs = 1367;                          % Solar Flux [W/m^2] *Varies based on distance to Sun*
gamma = pi/4;                       % Angle between orbit and satellite normal vector [rad]
                                    % Gamma is assumed to be between 0 and
                                    % pi/2 radians
inc = pi/180;                       % Inclination Angle [rad]
zeta = 0;
T_ambient = 302;                    % From Previous Space Team Calculations [K]
T_s = 4;                            % Temperature of Deep Space [K]
static = true;
PF = .2;                            % Reflectance of VED (percentage of flux)
EF = .2;                            % Reflectance of VED (percentage of flux)

% CONSTANTS
sigma = 5.6703 * 10^(-8);           % Stefan-Boltzmann Constant [W/m^2*K^4] 
rE = 6378140;                       % Earth radius [m]
r = rE + h;                         % Total Satellite Altitude [m]
R = r/rE;                           % Ratio of Radii
Q_earth = Cs * (1-albedo)/4;        % Planetary Infrared Flux [W/m^2]
Q_pla = Q_earth * (rE/r)^2;
D_flux = dissipation/A_tot;         % Dissipation Flux
M = 5.972 * 10^24;                  % Mass of Earth [kg]
G = 6.67 * 10^(-11);                % Gravitational Constant
w = 80 * 10^(-3);                   % Thickness of Plate [m]


% SHADOW AND PERIOD CALCULATIONS
P = round(2*pi * (r^3 / (G*M))^(.5));  % Orbital period of satellite (circular)
w = (2 * pi)/(P);                      % Satellite orbital angular velocity [rad/sec]
w_Sat = 2;                             % Satellite Rotational Angular Velocity [rad/sec]


% Variable Storage for Plotting
qVED = [];          % All Heat Fluxes into VED Side
qNadir = [];        % All Heat Fluxes into North Side
qNorth = [];        % All Heat Fluxes into Anti-Velocity Side
qSouth = [];        % All Heat Fluxes into Nadir Side
qVel = [];          % All Heat Fluxes into South Side
qAntiVel = [];      % All Heat Fluxes into Velocity Side
final_temp = [];

Qsun = [];
Qir = [];
Qalbedo = [];
zeta_Store = [];

T_VED = [T_ambient];
T_Nadir = [T_ambient];
T_North = [T_ambient];
T_South = [T_ambient];
T_Vel = [T_ambient];
T_AntiVel = [T_ambient];
time = -1:(2*P);

for t = 0:time(end)
    F_VED = sin(gamma)*cos(inc);
    if zeta >= 0 && zeta < (pi/2)
        F_North = 0;
        F_AntiVel = 0;
        F_Nadir = 0;
        F_South = cos(gamma)*cos(inc)*cos(zeta);
        F_Vel = cos(gamma)*cos(inc)*sin(zeta);  
    elseif zeta >= (pi/2) && zeta < pi
        F_AntiVel = 0;
        F_North = cos(gamma)*cos(inc)*sin(zeta);
        F_Nadir = 0;
        F_South = 0;
        F_Vel = cos(gamma)*cos(inc)*cos(zeta);  
    elseif zeta >= pi && zeta < (1.5*pi)
        F_South = 0;
        F_North = cos(gamma)*cos(inc)*cos(zeta);
        F_AntiVel = cos(gamma)*cos(inc)*sin(zeta);
        F_Nadir = 0;
        F_Vel = 0;  
    else
        F_South = cos(gamma)*cos(inc)*sin(zeta);
        F_North = 0;
        F_AntiVel = cos(gamma)*cos(inc)*cos(zeta);
        F_Nadir = 0;
        F_Vel = 0;  
    end
     
    % Solar Power Calculations
       
    Qsun_VED = PF * F_VED * Cs * absorptivity;
    Qsun_Nadir = EF * F_Nadir * Cs * absorptivity;
    Qsun_North = EF * F_North * Cs * absorptivity;
    Qsun_South = EF * F_South * Cs * absorptivity;
    Qsun_Vel = EF * F_Vel * Cs * absorptivity;
    Qsun_AntiVel = EF * F_AntiVel * Cs * absorptivity;
    
    Qsun = [Qsun; Qsun_VED, Qsun_Nadir, Qsun_North, Qsun_South, Qsun_Vel, Qsun_AntiVel];
   
    % Albedo View Factor Calculations
    zeta_North = angle_valid(pi + zeta);
    zeta_AntiVel = angle_valid(pi/2 + zeta);
    zeta_Vel = angle_valid(zeta - pi/2);
    zeta_South = angle_valid(zeta);
    
    zeta_Store = [zeta_Store; zeta_Vel, zeta_North, zeta_South, zeta_AntiVel];
    
    F_VED = 0;
    F_North = viewFactorAlbedo(zeta_North, h);
    F_AntiVel = viewFactorAlbedo(zeta_AntiVel, h);
    F_Nadir = 0;
    F_South = viewFactorAlbedo(zeta_South, h);
    F_Vel = viewFactorAlbedo(zeta_Vel, h);  
    
    %Qalbedo_VED = albedo * Cs *  F_VED;
    %Qalbedo_North = albedo * Cs * F_North;
    %Qalbedo_Vel = albedo * Cs * F_Vel;
    %Qalbedo_Nadir = albedo * Cs * F_Nadir;
    %Qalbedo_South = albedo * Cs * F_South;
    %Qalbedo_AntiVel = albedo * Cs * F_AntiVel;
    
    Qalbedo_VED = 0;
    Qalbedo_North = 0;
    Qalbedo_Vel = 0;
    Qalbedo_Nadir = 0;
    Qalbedo_South = 0;
    Qalbedo_AntiVel = 0;
    
    Qalbedo = [Qalbedo; Qalbedo_VED, Qalbedo_North, Qalbedo_South, Qalbedo_Nadir, Qalbedo_Vel, Qalbedo_AntiVel];

    F_Nadir = viewFactorAlbedo(gamma, h);
    % Planetary Infrared Flux Calculations
    Qir_VED = PF * Q_pla * emissivity * F_VED; 
    Qir_North = EF * Q_pla * emissivity * F_North; 
    Qir_South = EF * Q_pla * emissivity * F_South; 
    Qir_Nadir = EF * Q_pla * emissivity * F_Nadir;
    Qir_Vel = EF * Q_pla * emissivity * F_Vel; 
    Qir_AntiVel = EF * Q_pla * emissivity * F_AntiVel;  
    
    Qir = [Qir; Qir_VED, Qir_North, Qir_South, Qir_Nadir, Qir_Vel, Qir_AntiVel];
    
    %Qcond_VED = k * 
    % Energy Balance Equations
    Q_VED = Qsun_VED + Qir_VED + Qalbedo_VED + D_flux - emissivity * sigma * ((T_VED(end))^4 - T_ambient^4);
    Q_Nadir = Qsun_Nadir + Qir_Nadir + Qalbedo_Nadir + D_flux - emissivity * sigma * (T_Nadir(end)^4 - T_ambient^4);
    Q_North = Qsun_North + Qir_North + Qalbedo_North + D_flux - emissivity * sigma * (T_North(end)^4 - T_ambient^4);
    Q_South = Qsun_South + Qir_South + Qalbedo_South + D_flux - emissivity * sigma * (T_South(end)^4 - T_ambient^4);
    Q_Vel = Qsun_Vel + Qir_Vel + Qalbedo_Vel + D_flux - emissivity * sigma * (T_Vel(end)^4 - T_ambient^4);
    Q_AntiVel = Qsun_AntiVel + Qir_AntiVel + Qalbedo_AntiVel + D_flux - emissivity * sigma * (T_AntiVel(end)^4 - T_ambient^4);
    total = Q_VED + Q_Nadir + Q_North + Q_South + Q_Vel + Q_AntiVel;
      
    T_VED = [T_VED; transient(Q_VED, T_VED(end), A_s)];
    T_Nadir = [T_Nadir; transient(Q_Nadir, T_Nadir(end), A_s)];
    T_North = [T_North; transient(Q_North, T_North(end), A_l)];
    T_South = [T_South; transient(Q_South, T_South(end), A_l)];
    T_Vel = [T_Vel; transient(Q_Vel, T_Vel(end), A_l)];
    T_AntiVel = [T_AntiVel; transient(Q_AntiVel, T_AntiVel(end), A_l)];
    
    
    if ~static
      zeta = zeta + w_Sat*1;
      if zeta >= (2*pi)
          zeta = zeta - 2*pi;
      end
    end
    %if Q_VED > 1
    %    fprintf("t = %i, QVED = %4.4f\n", t, Q_VED);
    %end
end

plot(time, T_VED, time, T_Nadir, time, T_North, time, T_South, time, T_Vel, time, T_AntiVel)
fprintf("TEMPERATURE CALCULATIONS\n")
fprintf("Side            Equilibrium Temperature\n")
fprintf("VED Side             %4.4f K | %4.4f C\n", T_VED(end), T_VED(end) - 273);
fprintf("Nadir Side           %4.4f K | %4.4f C\n", T_Nadir(end), T_Nadir(end) - 273);
fprintf("North Side           %4.4f K | %4.4f C\n", T_North(end), T_North(end) - 273);
fprintf("South Side           %4.4f K | %4.4f C\n", T_South(end), T_South(end) - 273);
fprintf("Velocity Side        %4.4f K | %4.4f C\n", T_Vel(end), T_Vel(end) - 273);
fprintf("Anti-Velocity Side   %4.4f K | %4.4f C\n", T_AntiVel(end), T_AntiVel(end) - 273);



function Tf = transient(Q_ext, Ti, A)

delt = 1;
Cp = 887;
rho = 2700;             % Density of Aluminum [kg/m^3]
w = 80 * 10^(-3);       % Thickness of Plate [m]
m = rho * A * w;        % Mass of Plate [kg]
Tf = Ti + (Q_ext)*(delt/(m*Cp));

end

function Tf = steady_state(Q_ext, Ti, A)

delt = 1;
Cp = 887;
rho = 2700;             % Density of Aluminum [kg/m^3]
w = 80 * 10^(-3);       % Thickness of Plate [m]
m = rho * A * w;        % Mass of Plate [kg]
Tf = Ti + (A*Q_ext)*(delt/(m*Cp));

end

function F = viewFactorAlbedo(zeta, h)

rE = 6378140;
r = rE + h;
H = r/rE;
phi = asin(1/H);

W1 = .5 * asin(sqrt(H^2 - 1)/(H * sin(zeta)));
W2 = (1/(2 * H^2)) * (cos(zeta) * acos((-sqrt(H^2 - 1)*cot(zeta))))-(1/(2*H^2))*(sqrt(H^2-1)*sqrt(1-(H^2)*cos(zeta)^2));

if zeta <= (pi/2 - phi)
    F = cos(zeta)/(H^2);
elseif zeta > (pi/2 - phi) && zeta <= (pi/2 + phi)
    F = .5 - (2/pi) * (W1 - W2);
else
    F = 0;
end
end

function new_angle = angle_valid(angle)
if angle >= (2*pi)
    new_angle = angle - 2*pi;
elseif angle < 0
    new_angle = abs(angle);
else
    new_angle = angle;
end
end


