function rho2_search()
%% general
filename = 'inputHexane_N2.mat';

dat = readtable('PRSV_data.xlsx');
Na = 6.02214076e23; % 2019 SI units precisely
k  = 1.380649e-23;  % 2019 SI units precisely
R  = Na*k;

name = char(dat{[1,5],1})
Tc = dat{[1,5],2};
pc = dat{[1,5],3};
w  = dat{[1,5],4};  
kappa_1 = dat{[1,5],5};
kappa_2 = dat{[1,5],6};
kappa_3 = dat{[1,5],7};
%TL = dat{i,8};

T  = [310.93, 344.26, 377.55, 410.95, 444.25];

data1 = readtable('T31093.xlsx');
data2 = readtable('T34426.xlsx');
data3 = readtable('T37755.xlsx');
data4 = readtable('T41095.xlsx');
data5 = readtable('T44425.xlsx');

p1 = data1{:,1};
p2 = data2{:,1};
p3 = data3{:,1}; 
p4 = data4{:,1};
p5 = data5{:,1};
pdat = [p1;p2;p3;p4;p5]*6895; % Pa

lengths = [length(p1), length(p2), length(p3), length(p4), length(p5)];
xdat = [data1{:,2:3};data2{:,2:3};data3{:,2:3};data4{:,2:3};data5{:,2:3}];


figure(1)
hold on
suma = 1;

colorlist = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
markerlist = {'x','o','^','.','v'};
for i = 1:length(T)
    h(i,1:2) = plot(xdat(suma:(suma+lengths(i)-1),:),pdat(suma:(suma+lengths(i)-1))/1e6,markerlist{i},'MarkerSize',15);
    suma = suma+lengths(i);
    set(h(i,1:2),'LineWidth',2,'Color',colorlist{i})
end


%% n_HEXANE
% Tc(1) = 507.5;
% pc(1) = 3010000;	
% w(1) = 0.299;
% % a = 3.8271398089725;
% % b = 0.0001090644800;
% % M1=86.18e-3;%hexane
% vc(1) = 370e-6;

%% Nitrogen
% Tc(2) = 126.2;
% pc(2) = 3390000;
% vc(2) = 89.8e-6;
% w(2) = 0.039;

%% n_Decane--Methane
% Tc(1) = 617.7;
% pc(1) = 2110000;
% vc(1) = 617e-6;
% w(1) = 0.4923;
% 
% Tc(2) = 190.56;
% pc(2) = 4599000;
% vc(2) = 98.6e-6;
% w(2) = 0.0115;

%%
Tr = T./Tc;
kappa_0 = 0.378893 + 1.4897153*w - 0.17131848*w.^2 + 0.0196544*w.^3;

kappa   = kappa_0 + (kappa_1 + kappa_2.*(kappa_3 - Tr).*(1 - sqrt(Tr))).*(1 + sqrt(Tr)).*(0.7 - Tr);    %PRSV EoS

for i = 1:2
    for j = 1:length(T)
        if Tr(i,j) >= 0.7
            kappa(i,j) = kappa_0(i);
        end
    end
end

%kappa = (0.37464+1.54226*w-0.26992*w^2);       % PR EoS
alpha = (1+kappa.*(1-sqrt(Tr))).^2;
a = 0.457235*R^2*Tc.^2./pc.*alpha;
b = 0.077796*R*Tc./pc;

%     di = 0.305*Na^(-2/3);
%     K = di*a*b^(2/3)/R/T(i);
%rhoc = 1/vc;
idx = 0;
for i = 1:length(T)
    x0 = [ 815.1, 11.6, 311.9, 7554.1];
    for j = 1:lengths(i)
        %x0 = xdat(idx+j,:);
        p = pdat(idx+j);
        pars = [R,T(i),p,a(:,i)',b'];
        options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
            'MaxFunctionEvaluations',1000000000,'MaxIterations',10000000,'FunctionTolerance',1e-10,...
            'StepTolerance',1e-10);
        rho  = lsqnonlin(@(x) model(x,pars),x0,[0,0,0,0],[],options);
        rhoV = rho(1) + rho(2);
        rhoL = rho(3) + rho(4);
        res = model(rho,pars);
       
        xN2 = [rho(3)/rhoL, rho(1)/rhoV];
        output(idx+j,:) = [xN2, res];
        x0 = rho;

    end
        idx = idx + lengths(i);
        disp('progress: ')
        idx/sum(lengths)
end
%x0 = [ 815.1, 11.6, 311.9, 7554.1];
%x0 = [10, 6000, 5100, 1000]; % [rhoV1, rhoV2, rhoL1, rhoL2]

% options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxIterations',100000,'MaxFunctionEvaluations',100000)
% rho = fsolve(@(x) model(x,pars),x0,options)
save(filename,'output')
%load(filename,'output')

suma = 1;

for i = 1:length(T)
    g(i,1:2) = plot(output(suma:(suma+lengths(i)-1),1:2),pdat(suma:(suma+lengths(i)-1))/1e6);
    suma = suma+lengths(i);
    set(g(i,1:2),'LineWidth',2,'Color',colorlist{i})
end 

labelist = {'$T = 310.93 \mathrm{K}$','$T = 344.26 \mathrm{K}$','$T = 377.55 \mathrm{K}$','$T = 410.95 \mathrm{K}$','$T = 444.25 \mathrm{K}$'};
for i = 1:length(T)
    legend([h(i,1)],labelist{i},'interpreter','latex')
end
ax = gca;
%ax.YGrid = 'on';
set(gca,'fontsize',30)
set(gca,'TickLabelInterpreter','latex')
xlabel('$x_{\mathrm{N}_2}$','interpreter','latex')
ylabel('$p/\mathrm{MPa}$','interpreter','latex')

end

function res = model(rho,pars)
R = pars(1);
T = pars(2);
p = pars(3);
a = pars(4:5);
b = pars(6:7);

rhoV = rho(1) + rho(2);
rhoL = rho(3) + rho(4);

aV = a(1)*(rho(1)^2)/(rhoV^2) + a(2)*(rho(2)^2)/(rhoV^2) + 2*sqrt(a(1)*a(2))*(rho(1)*rho(2))/(rhoV^2);
aL = a(1)*(rho(3)^2)/(rhoL^2) + a(2)*(rho(4)^2)/(rhoL^2) + 2*sqrt(a(1)*a(2))*(rho(3)*rho(4))/(rhoL^2);

bV = b(1)*rho(1)/rhoV + b(2)*rho(2)/rhoV;
bL = b(1)*rho(3)/rhoL + b(2)*rho(4)/rhoL;

miV1 = R*T*log(abs(rho(1)/(1-bV*rhoV))) + R*T*b(1)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(1)*rho(1)+2*sqrt(a(1)*a(2))*rho(2))/(bV*rhoV) - aV*b(1)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(1)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miV2 = R*T*log(abs(rho(2)/(1-bV*rhoV))) + R*T*b(2)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(2)*rho(2)+2*sqrt(a(1)*a(2))*rho(1))/(bV*rhoV) - aV*b(2)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(2)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miL1 = R*T*log(abs(rho(3)/(1-bL*rhoL))) + R*T*b(1)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(1)*rho(3)+2*sqrt(a(1)*a(2))*rho(4))/(bL*rhoL) - aL*b(1)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(1)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));
miL2 = R*T*log(abs(rho(4)/(1-bL*rhoL))) + R*T*b(2)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(2)*rho(4)+2*sqrt(a(1)*a(2))*rho(3))/(bL*rhoL) - aL*b(2)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(2)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));

res(1) = miV1 - miL1;
res(2) = miV2 - miL2;

fV = R*T*rho(1)*(log(abs(rho(1)))-1) + R*T*rho(2)*(log(abs(rho(2)))-1) - R*T*rhoV*log(abs(1-bV*rhoV)) + (aV*rhoV)/2/sqrt(2)/bV*log(((1+(1-sqrt(2))*(bV*rhoV))/(1+(1+sqrt(2))*(bV*rhoV))));
fL = R*T*rho(3)*(log(abs(rho(3)))-1) + R*T*rho(4)*(log(abs(rho(4)))-1) - R*T*rhoL*log(abs(1-bL*rhoL)) + (aL*rhoL)/2/sqrt(2)/bL*log(((1+(1-sqrt(2))*(bL*rhoL))/(1+(1+sqrt(2))*(bL*rhoL))));

res(3) = p + fV - rho(1)*miV1 - rho(2)*miV2; %tlak plynu
res(4) = p + fL - rho(3)*miL1 - rho(4)*miL2; % tlak kvapaliny

end

function mi = potentials(rho,pars)
R = pars(1);
T = pars(2);
a = pars(3:4);
b = pars(5:6);

rhoV = rho(1) + rho(2);
rhoL = rho(3) + rho(4);

aV = a(1)*(rho(1)^2)/(rhoV^2) + a(2)*(rho(2)^2)/(rhoV^2) + 2*sqrt(a(1)*a(2))*(rho(1)*rho(2))/(rhoV^2);
aL = a(1)*(rho(3)^2)/(rhoL^2) + a(2)*(rho(4)^2)/(rhoL^2) + 2*sqrt(a(1)*a(2))*(rho(3)*rho(4))/(rhoL^2);

bV = b(1)*rho(1)/rhoV + b(2)*rho(2)/rhoV;
bL = b(1)*rho(3)/rhoL + b(2)*rho(4)/rhoL;
miV1 = R*T*log(abs(rho(1)/(1-bV*rhoV))) + R*T*b(1)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(1)*rho(1)+2*sqrt(a(1)*a(2))*rho(2))/(bV*rhoV) - aV*b(1)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(1)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miV2 = R*T*log(abs(rho(2)/(1-bV*rhoV))) + R*T*b(2)*rhoV/(1-bV*rhoV) + 1/2/sqrt(2)*((2*a(2)*rho(2)+2*sqrt(a(1)*a(2))*rho(1))/(bV*rhoV) - aV*b(2)/bV^2)*log(abs((1+(1-sqrt(2))*bV*rhoV)/((1+(1+sqrt(2))*bV*rhoV)))) - aV*rhoV*b(2)/(bV*(1-(rhoV*bV)^2+2*bV*rhoV));
miL1 = R*T*log(abs(rho(3)/(1-bL*rhoL))) + R*T*b(1)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(1)*rho(3)+2*sqrt(a(1)*a(2))*rho(4))/(bL*rhoL) - aL*b(1)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(1)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));
miL2 = R*T*log(abs(rho(4)/(1-bL*rhoL))) + R*T*b(2)*rhoL/(1-bL*rhoL) + 1/2/sqrt(2)*((2*a(2)*rho(4)+2*sqrt(a(1)*a(2))*rho(3))/(bL*rhoL) - aL*b(2)/bL^2)*log(abs((1+(1-sqrt(2))*bL*rhoL)/((1+(1+sqrt(2))*bL*rhoL)))) - aL*rhoL*b(2)/(bL*(1-(rhoL*bL)^2+2*bL*rhoL));

mi = [miV1, miV2, miL1, miL2];
end
