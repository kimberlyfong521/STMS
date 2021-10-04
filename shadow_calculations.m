clear
clc

d_E = 12756e3; %Diameter of the Earth (meters)
d_S = 2 * 695508e3; %Diameter of the Sun (meters)
d_ES = 149600000e3; %Distance from the Earth to the Sun (meters)

x_u = (d_ES * d_E)/(d_S - d_E);
x_p = (d_ES * d_E)/(d_S + d_E);

M = 5.972 * 10^24;              %mass of Earth (kg)
G = 6.67 * 10^(-11);            %Gravitational Constant

R = 300e3; %Altitude (meters)
R_total = R + d_E;

theta = acos(d_E/(2 * R_total)); 
alpha_u = asin(d_E/(2 * x_u));
alpha_p = asin(d_E/(2 * x_p));

dtheta = pi/1000;

beta = 0:dtheta:(2*pi);

period = round(2*pi * (R_total^3 / (G*M))^(.5));  %Orbital period of satellite (circular)
w = (2 * pi)/(period);           %Satellite angular velocity (rad/sec)

penumbra = [];
umbra = [];
no_shadow = [];
for i = beta
    if i >= (pi/2) && i <= (3 * pi/2)
        kappa = (x_p + R_total * cos(pi - i)) * tan(alpha_p);
        zeta  = (x_u - R_total * cos(pi - i)) * tan(alpha_u);
        
        delta = R_total * sin(abs(pi - i));
        rS = R_total * cos(abs(pi - i));
        if delta == kappa
            theta_k = beta;
        elseif delta == zeta
            theta_z = beta;
        end
        if delta <= kappa && delta > zeta
            shadow = 'Penumbra';
            %fprintf("Beta: %4.4f, Shadow: %s, Kappa = %4.2f, Zeta = %4.2f, Delta = %4.2f \n", i * 180/pi, shadow, kappa, zeta, delta);
            penumbra = [penumbra; i];
        elseif delta <= kappa && delta <= zeta
            shadow = 'Umbra';
            %fprintf("Beta: %4.4f, Shadow: %s, Kappa = %4.2f, Zeta = %4.2f, Delta = %4.2f \n", i * 180/pi, shadow, kappa, zeta, delta);
            umbra = [umbra; i];
        else
            shadow = 'No Shadow';
            %fprintf("Beta: %4.4f, Shadow: %s , Kappa = %4.2f, Zeta = %4.2f, Delta = %4.2f \n", i * 180/pi, shadow, kappa, zeta, delta);
            no_shadow = [no_shadow; i];
        end

    end
    
end

penumbra_entry = [];
penumbra_exit = [];
for j = 1:length(penumbra)
    if penumbra(j) <= umbra(j)
        penumbra_entry = [penumbra_entry, penumbra(j)];
    else
        penumbra_exit = [penumbra_exit, penumbra(j)];
        
    end
end

t_umbra = (umbra(end) - umbra(1))/w;
t_penEntry = (penumbra_entry(end) - penumbra_entry(1))/w;
t_penExit = (penumbra_exit(end) - penumbra_exit(1))/w;

fprintf("Time in Penumbra (Entry): %4.4f s\n", t_penEntry);
fprintf("Time in Umbra: %4.4f s\n", t_umbra);
fprintf("Time in Penumbra (Exit): %4.4f s\n", t_penExit);

