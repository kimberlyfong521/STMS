clear
clc
%% KEY ASSUMPTIONS
% 1. Earth is modeled as perfectly spherical
% 2. Satellite is CURRENTLY modeled to be in circular orbit
% 3. Polar Low-Earth Orbit
% 4. Orbit starts at the equator and follows polar path

%%
% CONSTANTS
albedo = 0.34;
emissivity = 0.85;
absorptivity = 0.13;
dissipation = 3;

altitude = 300 * 10^3;          %altitude (m)
M = 5.972 * 10^24;              %mass of Earth (kg)
G = 6.67 * 10^(-11);            %Gravitational Constant
r_earth = 6371 * 10^3;          %radius of Earth (m)
r = altitude + r_earth;         %Total orbital radius (m)

A0 = 0.03142;                                     % [m^2] Cross Sectional Area
As = A0*6;                                     % [m^2] Surface Area
sigma = 5.6703e-8;                             % [W/m^2*K^4] Stefan Constant
solar_flux = 1380;                             % [W/m^2] Solar Constant

% CALCULATIONS
period = round(2*pi * (r^3 / (G*M))^(.5));  %Orbital period of satellite (circular)
w = (2 * pi)/(period);           %Satellite angular velocity (rad/sec)
theta0 = pi/2 + acos(r_earth/r);            %Initial angle of eclipse
thetaf = 3*pi/2 - theta0;                %Final angle of eclipse  
temp_vals = [];

c = 4;
x = 1;

for t = 0:1:(c*period)
    zeta = w * t;       % Angle with respect to equatorial plane
    angle = zeta;
    if (zeta >= (x * theta0)) && (zeta <= (x * thetaf))
        solar = 0;
    else
        solar = 1;
    end
    if (zeta > (pi/2 + (x-1) * 2 * pi)) && (zeta <= (3 * pi/2 + (x-1) * 2 * pi))
        zeta = pi/2;
    end
    temp = avgtemp(albedo,emissivity,absorptivity,dissipation, solar, zeta, altitude);
    %fprintf("t = %i: zeta = %4.4f, T = %4.4f\n", t, zeta, temp);
    temp_vals = [temp_vals; t, temp, angle];
    
    if (t - period * x) == 0
        x = x + 1;
    end
end
sat_temp = temp_vals(1);




derivative = [];
for i = 1:length(temp_vals)
    if i < length(temp_vals)
        delta = temp_vals(i+1) - temp_vals(i);
        derivative = [derivative; delta];
    else
        derivative = [derivative; mean(derivative)];
    end
    
end

figure;
avg_flux = mean(derivative)/period;

fprintf("Avg Temperature Flux: %4.9f degrees K/sec\n", avg_flux);
fprintf("Avg Temperature Flux: %4.9f degrees K/hr\n", avg_flux*3600);
subplot(1, 3, 1)
plot(temp_vals(:, 1), temp_vals(:, 2), '-r')

xlabel("Time (sec)")
ylabel("Temperature (K)")

subplot(1, 3, 2)
plot(temp_vals(:, 3), temp_vals(:, 2), '-b')
xlabel("Angle (radians)")
ylabel("Temperature (K)")

subplot(1, 3, 3)
plot(temp_vals(:, 3) * 180/pi, temp_vals(:, 2), '-g')
xlabel("Angle (degrees)")
ylabel("Temperature (K)")


%plot(temp_vals(:, 1), derivative/3600, '-b')



