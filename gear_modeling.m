clc; clear all; close all

global Psi  Carrier_Time alpha_s alpha_r b_q delta_s delta_r Psi_sd Psi_rd Psi_sb Psi_rb delta_cr delta_cp N ks kr
global Omega_c Mass_sun Mass_ring Mass_carrier Mass_planet I_carrier I_sun I_ring I_planet
global D_Out_sun D_Out_planet D_Out_ring
global k_qp_sun k_qp_ring k_B_ring k_B_carrier k_B_planet k_cr T_m omega_c L

%% Define parameters in Example ( Table 1, Table 2, Table 3 of the paper)
N_T_sun = 16;                               % Number of Teeth of Sun
N_T_planet = 26;                            % Number of Teeth of Planet
N_T_ring = 68;                              % Number of Teeth of Ring

D_Out_sun = 254.5 / 1000;                   % Outer Diameter of Sun (m)
D_Out_planet = 383 / 1000;                  % Outer Diameter of Planet (m)
D_Out_ring = 980 / 1000;                    % Outer Diameter of Ring (m)

D_Root_sun = 202 / 1000;                    % Root Diameter of Sun (m)
D_Root_planet = 329 / 1000;                 % Root Diameter of Planet (m)
D_Root_ring = 930.6 / 1000;                 % Root Diameter of Ring (m)

D_Pitch_sun = 224 / 1000;                   % Pitch Diameter of Sun (m)
D_Pitch_planet = 364 / 1000;                % Pitch Diameter of Planet (m)
D_Pitch_ring = 952 / 1000;                  % Pitch Diameter of Ring (m)

T_Tooth_sun = 25 / 1000;                    % Transverse Tooth Thickness of Sun (m)
T_Tooth_planet = 18.5 / 1000;               % Transverse Tooth Thickness of Planet (m)
T_Tooth_ring = 25 / 1000;                   % Transverse Tooth Thickness of Ring (m)

Module = 14 / 1000;                         % Module of Sun, Planet and Ring (m)
Facewidth = 180 / 1000;                     % Facewidth of Sun, Planet and Ring (m)
Backlash = 0.482 / 1000;                    % Backlash of Sun-Planet and Planet-Ring (m)
Pres_angle = 24.6/180*pi;                   % Pressure Angle of Sun-Planet and Planet-Ring (m)
Gap = 0.526 / 1000;                         % Tooth Radial Gap of Sun-Planet and Planet-Ring (m)
Dist = 294 / 1000;                          % Sun/Pinion Center Distance (mm)

Mass_sun = 51;                              % Mass of Sun (kg)
Mass_ring = 28125;                          % Mass of Ring (kg)
Mass_carrier = 1330;                        % Mass of Carrier (kg)
Mass_planet = 114;                          % Mass of Planet (kg)

I_sun = 61.1;                               % Inertia of Sun (kg.m2)
I_ring = 2484;                              % Inertia of Ring (kg.m2)
I_carrier = 314.7;                          % Inertia of Carrier (kg.m2)
I_planet = 51.9;                            % Inertia of Planet (kg.m2)

k_qp_sun = 3.55*10^9;                       % refer Table 3
k_qp_ring = 4.56*10^9;                      % refer Table 3

k_B_ring = 126*10^6;                        % refer Table 3
k_B_carrier = 3.95*10^9;                    % refer Table 3
k_B_planet = 5.29*10^9;                     % refer Table 3
k_cr = 3.59*10^9;                           % refer Table 3

delta_cr = 1 / 1000;                        % refer Table 3 (m)
delta_cp = 0;                               % refer Table 3

k_u_sun = 3.02*10^6;                        % refer Table 3
k_u_ring = 24.4*10^6;                       % refer Table 3
k_u_carrier = 0;                            % refer Table 3
k_u_planet = 0;                             % refer Table 3

N = 3;                                      % Number of planets

%% Main part 
%% calc ks, kr ( reference parker, 2003, Mesh Phasing Relationships in Planetary and Epicyclic Gears )
% setting value
alpha1 = Pres_angle;
alpha2 = Pres_angle;
Rso = D_Out_sun / 2;
Rsb = D_Root_sun / 2;
Rpo = D_Out_planet / 2;
Rpb = D_Root_planet / 2;
Rro = D_Out_ring / 2;
Rrb = D_Root_ring / 2;
Zs = N_T_sun;
Zp = N_T_planet;
Zr = N_T_ring;
t_b = T_Tooth_planet;
omega_m = 34;
T_m = 2 * pi / omega_m;

Psi = linspace(0,2*pi,N+1);
Psi = Psi(1:N);
Gamma_s = Zs * Psi / (2*pi) - floor(Zs * Psi / (2*pi));                         % parker2003  (1)
Gamma_r = - (Zr * Psi / (2*pi) - floor(Zr * Psi / (2*pi)));                     % parker2003  (1)

B1E1 = sqrt(Rso^2 - Rsb^2) + sqrt(Rpo^2 - Rpb^2) - (Rsb + Rpb) * tan(alpha1);   % parker2003  (13)
B1P1 = Rsb * tan(alpha1) - (sqrt(Rso^2 - Rsb^2) - B1E1);                        % parker2003  (14)
p = 2*pi*Rsb/Zs;                                                                % parker2003  (15)
B1C1 = B1E1 - p;                                                                % parker2003  (15)
B1D1 = p;                                                                       % parker2003  (15)

O1O2 = (Rsb + Rpb) / cos(alpha1);                                               % parker2003  (16)
B2E2 = O1O2 * sin(alpha2) + sqrt(Rpo^2 - Rpb^2) - sqrt(Rro^2 - Rrb^2);          % parker2003  (16)
B2P2 = Rrb * tan(alpha2) - sqrt(Rro^2 - Rrb^2);                                 % parker2003  (17)
B2C2 = B2E2 - p;                                                                % parker2003  (18)
B2D2 = p;                                                                       % parker2003  (18)

Q2B2 = Rpb * tan(alpha1) + Rpb * (pi - alpha1 - alpha2) + Rpb * tan(alpha2)- B2P2 - t_b; % parker2003  (19)
B2Q3 = p * (1 - (Q2B2 / p - floor(Q2B2 / p)));                                           % parker2003  (20)
P2Q3 = abs(B2Q3 - B2P2);                                                                 % parker2003  (20)
gamma_rs = P2Q3 / p;                                                                     % parker2003  (21)

CT = 2;        % uint mesh cycle t / Tm
d_mt = 0.001;  % interval of mesh cycle
mesh_time = (0:d_mt:CT); % 
figure()

% initialize of ks and kr : the mesh tooth variation functions
ks = zeros(N, length(0:d_mt:CT));
kr = zeros(N, length(0:d_mt:CT));

% figure of Ksn, Krn
for j=1:N
    ks_y = (0:d_mt:CT);
    kr_y = (0:d_mt:CT);
    for i = 1: length(0:d_mt:CT)
        ks_y(i) = Ks(mesh_time(i), Gamma_s(j), p, B1D1, B1P1, B1C1);
        kr_y(i) = Ks(mesh_time(i), Gamma_r(j), p, B2D2, B2P2, B2C2);
    end
    ks(j, :) = ks_y;
    kr(j, :) = kr_y;
    
    % figure of Ksn (parker2003: Figure 2)
    subplot(N,2,j * 2 - 1)
    plt = plot(mesh_time, ks_y);
    plt(1).LineWidth = 2;
    xlabel('t/Tm mesh cycle')
    ylabel(['Ks' num2str(j,'%d') '-Teeth'])
    grid on
    
    % figure of Ksn (parker2003: Figure 2
    subplot(N,2,j * 2)
    plt = plot(mesh_time, kr_y);
    plt(1).LineWidth = 2;
    xlabel('t/Tm mesh cycle')
    ylabel(['Kr' num2str(j,'%d') '-Teeth'])
    grid on
end

ks = ks * k_qp_sun;
kr = kr * k_qp_ring;



%% Assumption about some parameters 
Psi = linspace(0,2*pi,N+1);
Psi = Psi(1:N);
alpha_s = Pres_angle;
alpha_r = Pres_angle;
b_q = 0.482;                           % backlashes - f. (7) 
delta_s = 0.5;
delta_r = 0.5;
Omega_c = omega_m / N_T_ring;
d_ct = d_mt;
Carrier_Time = (0:d_ct:CT);

Psi_sd = Psi - alpha_s;         % f. (1)
Psi_rd = Psi + alpha_r;         % f. (1)

Psi_sb = Psi + alpha_s;         % f. (5)
Psi_rb = Psi - alpha_r;         % f. (6)


zeta = 0.1;                     % Damping ratio
Omega = diag([42.202, 42.202, 69.491, 301.46, 305.59, 305.59, 700.06, 700.06, 1022.4, 1134.4, 1134.4, 1257.0, 1830.0, 1830.0, 1983.0, 2431.1, 2431.2, 2769.9]);                 % Natural frequency
omega_c = Omega(6,6);
L = 10^-6;

%% Define matrix U
M = zeros(9+3*N);
M(1:3,1:3) = diag([Mass_carrier, Mass_carrier, I_carrier/(D_Out_ring / 2)^2]);
M(4:6,4:6) = diag([Mass_ring, Mass_ring, I_ring/(D_Out_ring / 2)^2]);
M(7:9,7:9) = diag([Mass_sun, Mass_sun, I_sun/(D_Out_sun / 2)^2]);
for i = 1:N
    M((1:3) +9+3*(i-1),(1:3) +9+3*(i-1)) = diag([Mass_planet, Mass_planet, I_planet/(D_Out_planet / 2)^2]);
end

U = sqrt(M);
C = U'*(2*zeta*Omega)*U / omega_c;

%% Solving of (33),(34) by Newton method

% http://sites.science.oregonstate.edu/math/home/programs/undergrad/CalculusQuestStudyGuides/ode/second/so_num/so_num.html

Initial = [0,0,0,...
           0.1, -1.2, 0.55,...
           0.75, 0.6, 0.3,...
           0, -0.45, 0.33,...
           0.2, -0.25, 0.15,...
           -0.4, 0.35, -0.24];


N_c = length(Initial);
Sol = zeros(length(Initial), length(Carrier_Time));
Sol1 = zeros(length(Initial), length(Carrier_Time));
Sol(1:N_c,1) = Initial;             % Initial states
Sol1(1:N_c,1) = 15 / 1000;            % Initial velocity for each coordinate
inv_M = inv(M);

for j = 2:length(Carrier_Time)
  Sol(1:N_c,j) = Sol(1:N_c,j - 1) + d_ct* Sol(1:N_c, j-1);
  Sol1(1:N_c,j) = (Sol1(1:N_c,j - 1)' + d_ct * (-C * Sol1(1:N_c, j - 1) + RHS(j, Sol(1:N_c,j-1)'))' * inv_M)';
end

%% Plot solutions
figure()
subplot(2,2,1)
plot(Sol(1,:), Sol(2,:))
xlabel('x_c')
ylabel('y_c')
grid on


subplot(2,2,2)
plot(Sol(4,:), Sol(5,:))
xlabel('x_r')
ylabel('y_r')
grid on

subplot(2,2,3)
plot(Sol(7,:), Sol(8,:))
xlabel('x_s')
ylabel('y_s')
grid on

figure 
subplot(2,2,1)
plot(Carrier_Time, Sol(3,:))
xlabel('t')
ylabel('u_c')
grid on
       

subplot(2,2,2)
plot(Carrier_Time, Sol(6,:))
xlabel('t')
ylabel('u_r')
grid on

subplot(2,2,3)
plot(Carrier_Time, Sol(9,:))
xlabel('t')
ylabel('u_s')
grid on

%% Plot planets
figure 
subplot(2,2,1)
plot(Sol(10,:), Sol(11,:))
xlabel('\xi_1')
ylabel('\eta_1')
grid on
       

subplot(2,2,2)
plot(Sol(13,:), Sol(14,:))
xlabel('\xi_2')
ylabel('\eta_2')
grid on

subplot(2,2,3)
plot(Sol(16,:), Sol(17,:))
xlabel('\xi_3')
ylabel('\eta_3')
grid on

%% Function
%  parker2003  (9)
%  Ksn(t) and Krn(t) are periodic at the mesh period .
function res = Ks(t, tao, p, BD, BP, BC)
    for i = -5:5
        if t - tao >= (i * p - BP + BC) / p  && t - tao < (i * p + BD - BP) / p
            res = 1;
            return
        elseif t - tao >= (i * p + BD - BP) / p && t - tao < (i * p + BD - BP + BC) / p
            res = 2;
            return
        end
    end
end

%% Function
function res = RHS(t, X)
    global Psi  Carrier_Time alpha_s alpha_r b_q delta_s delta_r Psi_sd Psi_rd Psi_sb Psi_rb delta_cr delta_cp N ks kr
    global Omega_c Mass_sun Mass_ring Mass_carrier Mass_planet I_carrier I_sun I_ring I_planet
    global D_Out_sun D_Out_planet D_Out_ring
    global k_qp_sun k_qp_ring k_B_ring k_B_carrier k_B_planet k_cr T_m omega_c L

    %% X - vector of parameters (1 * (9 + 3N) vector) - (x_c,y_c,u_c,..., x_n, y_n, u_n)
    d = length(X);
    
    % Define additional parameters
    delta_sd = -X(7)*sin(Psi_sd) + X(8)*cos(Psi_sd)...
                -X(10:3:d)*sin(alpha_s) - X(11:3:d)*cos(alpha_s)...
                +X(12:3:d) + X(9);                                          % f.(1)
    delta_rd = -X(4)*sin(Psi_rd) + X(5)*cos(Psi_rd)...
                -X(10:3:d)*sin(alpha_r) - X(11:3:d)*cos(alpha_r)...
                -X(12:3:d) + X(6);                                          % f.(1)        
    
    %h_sd =   delta_sd>0;                                                   % f.(3)
    %h_rd =   delta_rd>0;                                                   % f.(3)
    
    %f_sd = k_qp_sun*ks(1,t)*h_sd.*delta_sd;                                % f.(2)
    %f_rd = k_qp_ring*kr(1,t)*h_rd.*delta_rd;                               % f.(2)
    
   
    delta_sb = -X(7)*sin(Psi_sb) - X(8)*cos(Psi_sb)...
                -X(10:3:d)*sin(alpha_s) + X(11:3:d)*cos(alpha_s)...
                -X(12:3:d) - X(9) ;                                         % f.(5)
    delta_rb =  X(4)*sin(Psi_rb) - X(5)*cos(Psi_rb)...
                +X(10:3:d)*sin(alpha_r) + X(11:3:d)*cos(alpha_r)...
                +X(12:3:d) - X(6);                                          % f.(6)        
              
    %h_sb =   delta_sb>b_q;                                                 % f.(8)
    %h_rb =   delta_rb>b_q;                                                 % f.(8)        
     
    %f_sb = k_qp_sun*ks(1, t)*h_sb.*(delta_sd - b_q);                       % f.(7)
    %f_rb = k_qp_ring*kr(1, t)*h_rb.*(delta_rd - b_q);                      % f.(7)
    
    delta_sR = X(7)*cos(Psi) + X(8)*sin(Psi) - X(10:3:d);                   % f. (9)
    delta_rR = -X(7)*cos(Psi) - X(8)*sin(Psi) + X(12:3:d);                  % f. (9)
    
    delta_sb_mod = delta_sb - 2*delta_s*sin(alpha_s);                       % f.(10)
    delta_rb_mod = delta_rb - 2*delta_r*sin(alpha_s);                       % f.(10)
    
    h_sw = delta_sR>delta_s;                                                % f.(12)
    h_rw = delta_rR>delta_r;                                                % f.(12)
    
    
    f_sd = k_qp_sun*ks(1, t)*h_sw.*delta_sd;                                % f. (11)
    f_rd = k_qp_ring*kr(1, t)*h_rw.*delta_rd;                               % f. (11)
    
    f_sb = k_qp_sun*ks(1, t)*h_sw.*delta_sb_mod;                            % f. (11)
    f_rb = k_qp_ring*kr(1, t)*h_rw.*delta_rb_mod;                           % f. (11)
    
    
    delta_c = ((X(1)*cos(Psi) + X(2)*sin(Psi)-X(10:3:d)).^2+...
        (-X(1)*sin(Psi) + X(2)*sin(Psi) + X(3)-X(11:3:d)).^2).^0.5;         %f. (13)
    
    theta_c = atan((-X(1)*sin(Psi) + X(2)*sin(Psi) + X(3)-X(11:3:d))./...
        (X(1)*cos(Psi) + X(2)*sin(Psi)- X(10:3:d)));                        % f. (14)

    mu_c = delta_c>delta_cp;                                                % f. (16)
    
    f_c_x =  0;
    f_c_y =  0;
    f_c_u =  0;

    f_cx = k_B_planet*(delta_c-delta_cp).*mu_c.*cos(theta_c+Psi);           % f. (15)
    f_cy = k_B_planet*(delta_c-delta_cp).*mu_c.*sin(theta_c+Psi);           % f. (15)
    f_cu = k_B_planet*(delta_c-delta_cp).*mu_c.*sin(theta_c);               % f. (15)
    f_pxi = -k_B_planet*(delta_c-delta_cp).*mu_c.*cos(theta_c);             % f. (15)
    f_peta = -k_B_planet*(delta_c-delta_cp).*mu_c.*sin(theta_c);            % f. (15)
     
    delat_cr_r = sqrt((X(1) - X(4))^2 + (X(2) - X(5))^2);                   % f. (17)
    delat_sr_r = sqrt((X(7) - X(4))^2 + (X(8) - X(5))^2);                   % f. (17)
    
    theta_cr = atan((X(2) - X(5))/(X(1) - X(4)));                           % f. (18)
    theta_sr = atan((X(8) - X(5))/(X(7) - X(4)));                           % f. (18)
    
    mu_cr = delat_cr_r > delta_cr;                                          % f. (20)
    mu_sr = delat_sr_r > delta_cr;                                          % f. (20)
    
    %k_B_ring k_B_carrier k_B_planet k_cr 
    f_crx = k_cr*(delat_cr_r - delta_cr)*mu_cr*cos(theta_cr);               % f. (19)
    f_cry = k_cr*(delat_cr_r - delta_cr)*mu_cr*sin(theta_cr);               % f. (19)
    
    f_srx = k_cr*(delat_sr_r - delta_cr)*mu_sr*cos(theta_sr);               % f. (19)
    f_sry = k_cr*(delat_sr_r - delta_cr)*mu_sr*sin(theta_sr);               % f. (19)
    
    res = zeros(9+3*N,1);
    res(4:6) = [sum(f_rd.*sin(Psi_rd)), -sum(f_rd.*cos(Psi_rd)), -sum(f_rd)];   % f. (27)
    res(7:9) = [sum(f_sd.*sin(Psi_sd)), -sum(f_sd.*cos(Psi_sd)), -sum(f_sd)];   % f. (27)
    
    for j = 1:N
        res((7+3*j):(9+3*j)) = [f_sd(j)*sin(alpha_s)-f_rd(j)*sin(alpha_r),...
                                f_sd(j)*cos(alpha_s)+f_sd(j)*cos(alpha_r),...
                                -f_sd(j)+f_rd(j)];                              %. f(27)
    end
    
    res(4:6) = res(4:6) + [-sum(f_rb.*sin(Psi_rb)), sum(f_rb.*cos(Psi_rb)), sum(f_rb)]';    %. f(28)
    res(7:9) = res(7:9) + [-sum(f_sb.*sin(Psi_sb)), sum(f_sb.*cos(Psi_sb)), sum(f_sd)]';    %. f(28)
    
    for j = 1:N
        res((7+3*j):(9+3*j)) =  res((7+3*j):(9+3*j)) + ...
            [f_sb(j)*sin(alpha_s)-f_rb(j)*sin(alpha_r),...
             -f_sb(j)*cos(alpha_s)-f_sb(j)*cos(alpha_r),...
             f_sb(j)-f_rb(j)]';                                                             %. f(28)
    end
    
    res(1:3) = [sum(f_cx)+f_crx + f_c_x, sum(f_cy)+f_cry + f_c_y, sum(f_cu) + f_c_u ]';     % f. (29) SOME PARAMETERS ARE UNDEFINED!!!!
    res(4:6) = res(4:6) + [f_crx + f_srx, f_cry + f_sry, 0]';                               % f. (29) SOME PARAMETERS ARE UNDEFINED!!!!
    res(7:9) = res(7:9) + [f_srx, f_sry,0]';                                                % f. (29) SOME PARAMETERS ARE UNDEFINED!!!!
    
    for j = 1:N
        % f. (29) SOME PARAMETERS ARE UNDEFINED!!!!
        res((7+3*j):(9+3*j)) =  res((7+3*j):(9+3*j)) + ...
            [f_pxi(j), f_peta(j), 0]';
    end
    
    % Define right hand side of equation 
    F  = zeros(9+3*N,1);
    
    % Omega_c Mass_sun Mass_ring Mass_carrier Mass_planet

    f_cgx = - Mass_carrier*sin(Omega_c*Carrier_Time(t)/Omega_c);            % f. (25)
    f_cgy = - Mass_carrier*cos(Omega_c*Carrier_Time(t)/Omega_c);            % f. (25)
    
    f_sgx = - Mass_sun*sin(Omega_c*Carrier_Time(t)/Omega_c);                % f. (25)
    f_sgy = - Mass_sun*cos(Omega_c*Carrier_Time(t)/Omega_c);                % f. (25)
    
    f_rgx = - Mass_ring*sin(Omega_c*Carrier_Time(t)/Omega_c);               % f. (25)
    f_rgy = - Mass_ring*cos(Omega_c*Carrier_Time(t)/Omega_c);               % f. (25)
    
    f_xigx = - Mass_planet*sin(Omega_c*Carrier_Time(t)/Omega_c + Psi);      % f. (26)
    f_etagy = - Mass_planet*cos(Omega_c*Carrier_Time(t)/Omega_c + Psi);     % f. (26)
    
    
    % F 
    F(1:3) = [f_cgx, f_cgy, 0]';
    F(4:6) = [f_rgx, f_rgy, 0]';
    F(7:9) = [f_sgx, f_sgy, 0]';
    
    F(10:3:d) = f_xigx;
    F(11:3:d) = f_etagy;
    
    res(1:3) = F(1:3) - res(1:3)./[Mass_carrier, Mass_carrier, I_carrier/(D_Out_ring/2)^2]';
    res(4:6) = F(4:6) - res(4:6)./[Mass_ring, Mass_ring, I_ring/(D_Out_ring/2)^2]';
    res(7:9) = F(7:9) - res(7:9)./[Mass_sun, Mass_sun, I_sun/(D_Out_sun/2)^2]';
    
    for j = 1:N
        res((7+3*j):(9+3*j)) = F((7+3*j):(9+3*j))...
            - res((7+3*j):(9+3*j))./[Mass_planet, Mass_planet, I_planet/(D_Out_planet/2)^2]';
    end
    
    res = res / (L * omega_c ^ 2);
end