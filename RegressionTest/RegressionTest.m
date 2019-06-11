function Passed=RegressionTest

% this is a quick synthetic test case comparing PICO-Ua with the analytical
% solutions for a perfectly circular ice shelf. Run by typing
% 'RegressionTest', should result in a Pass!

clearvars;


%% Run Ua on a quadratic domain with a round ice shelf at its center
% Adjust the number of elements in Ua2D_InitialUserInput to increase 
% accuracy / speed
Ua

%% Run PICO 
% adjust further options in RunPICO

CallToPICO('watershed'); % 'polygon';

%% Compare to direct solution

% load Ua Output
load ./ResultsFiles/0000001-RegressionTest.mat;
% load PICO results
load ./ResultsFiles/0000001-RegressionTest_PICO.mat


% Geometric specifications 
number_of_boxes = PICO_opts.nmax;

% Use a theoretically derived area of boxes (see Appendix B in Paper)
radius = 500*1e3; %m
circle_area = pi*radius^2; % m^2 
area_per_box = circle_area/number_of_boxes; %m^2

ice_thickness = median(F.h); % m
rhoi= median(F.rho); % kg m^3
rhow = F.rhow; % kg m^3

% parameters
gammaT = 2e-5;
C = 1e6;

a1 = -0.0572; % deg C PSU-1
b1 = 0.0788; % deg C
c1 = 7.77e-8; % deg C Pa-1
alph = 7.5e-5; % deg C-1
beta = 7.7e-4; % PSU-1
rhostar = 1033;
L = 3.34e5; %J Kg-1
cp = 3974; %J Kg-1 degC-1

nu = 900/1000; %rho_i / rhoo
lambda = L/cp;
g = 9.81;

Tkm_test = zeros(size(Tkm));
Skm_test = zeros(size(Skm));
Mk_test = zeros(size(Tkm)); %uniform in each box
q_test = 0; % uniform per box

% solution box 1

% T1, S1, m1, q
g1 = area_per_box * gammaT;
g2 = g1/ (nu*lambda);
s = S0 / (nu*lambda);

pressure = rhoi*g*ice_thickness;
Tstar = a1*S0 + b1 - c1*pressure - T0; 

tmp = g1/(2*C *rhostar * (beta*s - alph));
x = - tmp + sqrt( tmp^2 - g1*Tstar/(C *rhostar * (beta*s - alph))  );
Tkm_test(1) = T0 - x;
y =  S0*x/(nu*lambda);
Skm_test(1) = S0 - y;
q_test = C*rhostar*(beta*(S0-Skm(1)) - alph*(T0-Tkm(1)));
Mk_test(1) = (- gammaT * (a1*Skm_test(1) + b1 - c1*pressure - Tkm_test(1))/(nu*lambda)) * 86400 * 365.25; %m per a


% solution other boxes
for k = 2:number_of_boxes
    
    g1 = area_per_box * gammaT;
    g2 = g1/ (nu*lambda);
    
    Tstar = a1*Skm_test(k-1) + b1 - c1*pressure - Tkm_test(k-1); 
    
    x = -g1*Tstar / (q_test + g1 - g2 *a1 * Skm_test(k-1));
    Tkm_test(k) = Tkm_test(k-1) - x; 

    y = Skm_test(k-1)*x / (nu *lambda);
    Skm_test(k) = Skm_test(k-1) - y;
    
    Mk_test(k) = - gammaT * (a1*Skm_test(k) + b1 - c1*pressure - Tkm_test(k))/(nu*lambda) * 86400 * 365.25;
end

%

Passed = 0;
threshold = 0.1;

if (max(abs(Mk_test - [median(Mk(PBOX==1)), median(Mk(PBOX==2)), median(Mk(PBOX==3)), median(Mk(PBOX==4)), median(Mk(PBOX==5)) ])) < threshold) & (abs(q-q_test) < threshold)
    
    Passed = 1;
end

fprintf(strcat('Basal melt difference: ', num2str(Mk_test - [median(Mk(PBOX==1)), median(Mk(PBOX==2)), median(Mk(PBOX==3)), median(Mk(PBOX==4)), median(Mk(PBOX==5)) ]), ' m/a \n'));
fprintf(strcat('overturning difference: ', num2str(q-q_test), ' \n'));
fprintf(strcat('Passed = ',num2str(Passed), '\n'));

end
