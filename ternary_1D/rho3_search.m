%function rho3_search()
close all, clear
%% general
Na=6.02214076e23;
k=1.380649e-23;
R=Na*k;
%T=320.15;
T=298.15;
p = 101325;
x3vec = linspace(0.1,0.7,20);
A = zeros(length(x3vec),3);
B = zeros(length(x3vec),3);
SSE = zeros(length(x3vec),1);
%% Water
Tc(1) = 647.286;    % Trange 274–623
pc(1) = 22089750;
w(1)  = 0.3438;
kappa_1(1) = -0.06635;
kappa_2(1) = 0.0199;
kappa_3(1) = 0.443;
%% Dodecane
Tc(2) = 658.2;      %Trange 312-520
pc(2) = 1823830;
w(2)  = 0.57508;
kappa_1(2) = 0.05426;
kappa_2(2) = 0.8744;
kappa_3(2) = 0.505;
%% 1-Butanol
Tc(3) = 562.98;     %Trange 352–399
pc(3) = 4412660;
w(3)  = 0.59022;
kappa_1(3) = 0.33431;
kappa_2(3) = -1.1748;
kappa_3(3) = 0.642;

%%
Tr = T./Tc;
kappa_0 = 0.378893 + 1.4897153*w - 0.17131848*w.^2 + 0.0196544*w.^3;

kappa  = kappa_0 + (kappa_1 + kappa_2.*(kappa_3 - Tr).*(1 - sqrt(Tr))).*(1 + sqrt(Tr)).*(0.7 - Tr);    %PRSV EoS
            
%kappa = (0.37464+1.54226*w-0.26992*w^2);       % PR EoS
alpha = (1+kappa.*(1-sqrt(Tr))).^2;
a = 0.457235*R^2*Tc.^2./pc.*alpha;
b = 0.077796*R*Tc./pc;
x0(1:3) = 1e-1;
x0(4:6) = 1e1;
x0 = x0*1e2;
x0 = [4e3, 4e3, 4e3, 1e4, 1e2, 1e2]; % [rhoV1, rhoV2, rhoV3, rhoL1, rhoL2, rhoL3]

for i = 1:length(x3vec)

    pars = [R,T,p,a,b,x3vec(i)];
    %x0 = [1.3044684, 43.878, 102.1195, 55345.5454, 4379.73346, 10846.982];

    residualy = model(x0,pars);
    residualy';

    %options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
    %rho = fsolve(@(x) model(x,pars),x0,options); 
    %rho'

    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
            'MaxFunctionEvaluations',1000000000,'MaxIterations',10000000,'FunctionTolerance',1e-10,...
            'StepTolerance',1e-10);
%     options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
%         'MaxFunctionEvaluations',1000000000,'MaxIterations',10000000,'FunctionTolerance',1e-10,...
%         'StepTolerance',1e-10);
    rho  = lsqnonlin(@(x) model(x,pars),x0,[0,0,0,0,0,0],[],options);
     
     residualy = model(rho,pars);
     SS = sum(residualy.^2);
     
     rho';
     rhoV = rho(1) + rho(2) + rho(3);
     rhoL = rho(4) + rho(5) + rho(6);
     y = rho(4:6)/rhoL;
     x = rho(1:3)/rhoV;
        
     A(i,:) = x;
     B(i,:) = y;
     SSE(i) = SS;
     mess = fprintf('Iteration %i done !',i);
     disp(mess)

end


disp(SSE)

figure
%-- Plot the axis system
[h,hg,htick]=terplot;
%-- Plot the data ...
hter=ternaryc(A(:,1),A(:,2),A(:,3));
%-- ... and modify the symbol:
set(hter,'marker','o','markerfacecolor','none','markersize',15,'color','#A2142F')
hlabels=terlabel('water','dodecane','butanol');
%end

function res = model(rho,pars)
R = pars(1);
T = pars(2);
p = pars(3);
a = pars(4:6);
b = pars(7:9);
x3fix = pars(10);

rhoV = rho(1) + rho(2) + rho(3);
rhoL = rho(4) + rho(5) + rho(6);

x = rho(4:6)/rhoL; % liquid mol. fractions
y = rho(1:3)/rhoV; % vapour mol. fractions

aV = a(1)*y(1)^2 + a(2)*y(2)^2 + a(3)*y(3)^2 + ...
      2*sqrt(a(1)*a(2))*y(1)*y(2) + ...
      2*sqrt(a(3)*a(2))*y(3)*y(2) ...
      + 2*sqrt(a(1)*a(3))*y(1)*y(3);
aL = a(1)*x(1)^2 + a(2)*x(2)^2 + a(3)*x(3)^2 + ...
      2*sqrt(a(1)*a(2))*x(1)*x(2) + ...
      2*sqrt(a(1)*a(3))*x(1)*x(3) ...
      + 2*sqrt(a(2)*a(3))*x(2)*x(3);

bV = b(1)*y(1) + b(2)*y(2) + b(3)*y(3);
bL = b(1)*x(1) + b(2)*x(2) + b(3)*x(3);

miV1 = R*T*log(abs(rho(1)/(1-bV*rhoV))) + R*T*b(1)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(1)*rho(1) + 2*sqrt(a(1)*a(2))*rho(2) + 2*sqrt(a(1)*a(3))*rho(3))/(bV*rhoV) - aV*b(1)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(1)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miV2 = R*T*log(abs(rho(2)/(1-bV*rhoV))) + R*T*b(2)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(2)*rho(2) + 2*sqrt(a(1)*a(2))*rho(1) + 2*sqrt(a(3)*a(2))*rho(3))/(bV*rhoV) - aV*b(2)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(2)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miV3 = R*T*log(abs(rho(3)/(1-bV*rhoV))) + R*T*b(3)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(3)*rho(3) + 2*sqrt(a(1)*a(3))*rho(1) + 2*sqrt(a(3)*a(2))*rho(2))/(bV*rhoV) - aV*b(3)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(3)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miL1 = R*T*log(abs(rho(4)/(1-bL*rhoL))) + R*T*b(1)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(1)*rho(4) + 2*sqrt(a(1)*a(2))*rho(5) + 2*sqrt(a(1)*a(3))*rho(6))/(bL*rhoL) - aL*b(1)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(1)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));
miL2 = R*T*log(abs(rho(5)/(1-bL*rhoL))) + R*T*b(2)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(2)*rho(5) + 2*sqrt(a(1)*a(2))*rho(4) + 2*sqrt(a(3)*a(2))*rho(6))/(bL*rhoL) - aL*b(2)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(2)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));
miL3 = R*T*log(abs(rho(6)/(1-bL*rhoL))) + R*T*b(3)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(3)*rho(6) + 2*sqrt(a(1)*a(3))*rho(4) + 2*sqrt(a(3)*a(2))*rho(5))/(bL*rhoL) - aL*b(3)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(3)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));
res(1) = miV1 - miL1;
res(2) = miV2 - miL2;
res(3) = miV3 - miL3;

fV = R*T*rho(1)*(log(abs(rho(1)))-1) + R*T*rho(2)*(log(abs(rho(2)))-1) + R*T*rho(3)*(log(abs(rho(3)))-1) - R*T*rhoV*log(abs(1-bV*rhoV)) + (aV*rhoV)/2/sqrt(2)/bV*log(((1+(1-sqrt(2))*(bV*rhoV))/(1+(1+sqrt(2))*(bV*rhoV))));
fL = R*T*rho(4)*(log(abs(rho(4)))-1) + R*T*rho(5)*(log(abs(rho(5)))-1) + R*T*rho(6)*(log(abs(rho(6)))-1) - R*T*rhoL*log(abs(1-bL*rhoL)) + (aL*rhoL)/2/sqrt(2)/bL*log(((1+(1-sqrt(2))*(bL*rhoL))/(1+(1+sqrt(2))*(bL*rhoL))));

res(4) = p + fV - rho(1)*miV1 - rho(2)*miV2 - rho(3)*miV3; %tlak plynu
res(5) = p + fL - rho(4)*miL1 - rho(5)*miL2 - rho(6)*miL3; % tlak kvapaliny
% P1 = fV - rho(1)*miV1 - rho(2)*miV2 - rho(3)*miV3; %tlak plynu
% P2 = + fL - rho(4)*miL1 - rho(5)*miL2 - rho(6)*miL3; % tlak kvapaliny
% res(4) = P1 - P2;
% res(5) = 1e3*(x(2) - 0.35);
res(6) = 1e3*(x(3) - x3fix);
%res(6) = rho(1)/rhoV + rho(2)/rhoV + rho(3)/rhoV - rho(4)/rhoL - rho(5)/rhoL - rho(6)/rhoL; 
end