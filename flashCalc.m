%% INITIALIZE MATLAB
clear all;
clc;
close all;
format long;

%% DASHBOARD
Z = [0.3805 0.0933 0.0885 0.0600 0.0378 .0356 .3043]; % Component Mole Fraction
P = 6.89476; % MPa
T = 327.594; % Kelvin
R = 8.31446261815324e-3; % L.MPa/K/mol 
MW = [16.043 30.070 44.096 58.123 72.150 86.177 200]; % kg/kgmol

% % Calculating C7+ Pseudocritical Values % %
gamma_c7 = .8366;
MW_c7 = 200;
TB = (4.5579*(MW_c7^.15178)*gamma_c7^.15427)^3; % Rankine
Ppc = exp(8.3634 - .0566/gamma_c7 ...
        - (.24244+2.2898/gamma_c7+.11857/(gamma_c7)^2)*1e-3*TB ...
        + (1.4685 + 3.648/gamma_c7 + .47227/(gamma_c7)^2)*1e-7*TB^2 ...
        - (.42019 + 1.6977/(gamma_c7)^2)*1e-10*TB^3); % psia
Tpc = 341.7 + 811*gamma_c7 + (.4244 + .1174*gamma_c7)*TB ...
    + (.4669 - 3.2623*gamma_c7)*1e5/TB; % Rankine
TBr = TB/Tpc;
Ppc_MPa = Ppc/145.04; % MPa
Tpc_Kelvin = (Tpc - 491.67)/1.8 + 273.15; % Kelvin
% ref: McCaine Eq. B-12 to B-14 

% % Calculating C7+ acentric factor % %
A1 = -5.92714; A2 = 6.09648; A3 = 1.28862; A4 = -0.169347; A5 = 15.2518;
A6 = -15.6875; A7 = -13.4721; A8 = 0.43577;
omega_c7 = (-log(Ppc/14.7) + A1 + A2/TBr + A3*log(TBr) + A4*TBr^6)/...
    (A5 + A6/TBr + A7*log(TBr) + A8*TBr^6);
% ref: Eq. 5.60 Phase Behavior Monograph Ch5, because TBr<0.8

% Critical Parameters and Acentric Factor
Pc = [4.599 4.872 4.248 3.796 3.370 3.025 Ppc_MPa]; % MPa
Tc = [190.56 305.32 369.83 425.12 469.7 507.6 Tpc_Kelvin]; % Kelvin
Pr = P./Pc;
Tr = T./Tc;
omega = [0.0115 0.0995 0.1523 0.2002 0.2515 0.3013 omega_c7];

%% 1. Calculating Ki
K = (1./Pr).*exp(5.37.*(1+omega).*(1-(1./Tr)));

% Checking to see if the fluid is in two-phase region.
if sum(K.*Z) < 1 
    fprintf('The mixture is a compressed liquid.');
    break
elseif sum(Z./K) < 1
    fprintf('The mixture is an undersaturated vapour.');
    break
end  
%% 2. Calculate x & y using Eq. 5.4-6
err1 = inf;
iter = 0;
while err1 > 1e-12  
% calculating beta by solving the Rachford-Rice equation; ref: Pederson eq.
% 6.22
beta = 0.5; % guess
err = inf;
while err>1e-6
F = sum(Z.*(K-1)./(1 + beta.*(K-1)));
dF = -sum((Z.*(K-1).^2)./((1+beta.*(K-1)).^2));
beta_old = beta;
beta = beta - F/dF;
err = abs(beta-beta_old);
end

x = Z./(1+(K-1).*beta);
y = K.*Z./(1+(K-1).*beta);

%% 3. PR EOS For Vapor phase
%%%%%%%%%%%% Vapor %%%%%%%%%%%%%
m = .375 + 1.574.*omega - .267.*omega.^2;
a_vector = .45724.*(R.^2).*(Tc.^2)./Pc.*(1+ m.*(1-sqrt(Tr))).^2;
b_vector = .07780*R*Tc./Pc;

% % Calculating a & b using Mixing Rule % %
% Binary Interaction Parameters; ref: Table A.4.3, Danesh (1998)
BIP     =         [  0 0 0 0 0 0 0;
                 .0026 0 0 0 0 0 0;
             .0140 .0011 0 0 0 0 0;
         .0133 .0096 .0033 0 0 0 0;
     .0236 .0078 .0120 .0170 0 0 0;
 .0422 .0140 .0110 .0240 .0174 0 0;
            .05 .03 .02 .001 0 0 0];
BIP = BIP + BIP';
        
b = sum(y.*b_vector);     
a = 0;        
Sum = 0;
for i = 1:length(Z)
    for j = 1:length(Z)
        Sum = y(i)*y(j)*((a_vector(i)*a_vector(j))^0.5)*(1-BIP(i,j));
        a = Sum + a;
    end
end
 
% % Calculating gas molar volume to obtain its z factor
syms vg
eqn_g = -P*vg^3 + (R*T-P*b)*vg^2 ...
    + (2*b*R*T - a + 3*P*b^2)*vg ...
    - P*b^3 - R*T*b^2 + a*b == 0;
S = solve(eqn_g, vg, 'Real', true);

% EQ = [-P, R*T-P*b, 2*b*R*T - a + 3*P*b^2, -P*b^3 -R*T*b^2 +a*b  ];

Vg = double(S); % Ltrs/mole
Zg = P*Vg/R/T;

%% 4. Calculating component Fugacity in gas phase 
A = a.*P./((R.*T).^2);
B = b.*P./(R.*T);

phi_G = 0;
for i = 1:7
    sigma = 0; Bi = b_vector(i)*P/R/T;
    for j = 1:7
        sigma =  sigma + y(j)*(1-BIP(i,j))*(a_vector(i)*a_vector(j))^0.5;
    end
    phi_G(i) = exp(Bi./B*(Zg-1)-log(Zg-B)+A./(B.*(-2*sqrt(2))).* ...
        (Bi./B-2/a*sigma)*log((Zg+(1-sqrt(2)).*B)./(Zg+(1+sqrt(2)).*B)));
end

fg = phi_G.*y.*P;

%% 5. PR EOS For Liquid phase
%%%%%%%%%%%%%% Liquid %%%%%%%%%%%%%

bl = sum(x.*b_vector);

% a for liquid phase
al = 0;        
Sum = 0;
for i = 1:length(Z)
    for j = 1:length(Z)
        Sum = x(i)*x(j)*((a_vector(i)*a_vector(j))^0.5)*(1-BIP(i,j));
        al = Sum + al;
    end
end
 
% Calculating liquid molar volume to obtain its z factor
syms vl
eqn_l = -P*vl^3 + (R*T-P*bl)*vl^2 ...
    + (2*bl*R*T - al + 3*P*bl^2)*vl ...
    - P*bl^3 - R*T*bl^2 + al*bl == 0;
Sl = solve(eqn_l, vl, 'Real', true);
Vl = double(Sl); % Ltrs/mole
Zl = P*Vl/R/T;

%% 6. Calculating component Fugacity in liquid phase 
A = al.*P./((R.*T).^2);
B = bl.*P./(R.*T);

phi_L = 0;
for i = 1:7
    sigma = 0; Bi = b_vector(i)*P/R/T;
    for j = 1:7
        sigma =  sigma + x(j)*(1-BIP(i,j))*(a_vector(i)*a_vector(j))^0.5;
    end
    phi_L(i) = exp(Bi./B*(Zl-1)-log(Zl-B)+A./(B.*(-2*sqrt(2))).* ...
        (Bi./B-2/al*sigma)*log((Zl+(1-sqrt(2)).*B)./(Zl+(1+sqrt(2)).*B)));
end

fl = phi_L.*x.*P;

%% 7. Determining Error Criteria and Calculating the New K
err1 = sum((1 - fl./fg).^2);
K = K.*(fl./fg);
iter = iter + 1;
end

%% Table of Results
sprintf(' The results after %d iterations are: ', iter)
disp(['gas phase molar amount (nV): ', num2str(beta)]);
disp(['liquid phase molar amount (nL): ', num2str(1-beta)]);

VarNames = {'GasComposition', 'LiquidComposition'};
names = {'C1', 'C2', 'C3', 'nC4', 'nC5', 'nC6', 'C7+'};
results = table(y', x');
results.Properties.RowNames = names;
results.Properties.VariableNames = VarNames;
disp(results);

% Reporting Molecular Weights
gasMW = sum(y.*MW);
liqMW = sum(x.*MW);

disp(['gas phase Molecular Weight: ', num2str(gasMW)]);
disp(['liquid phase Molecular Weight: ', num2str(liqMW)]);