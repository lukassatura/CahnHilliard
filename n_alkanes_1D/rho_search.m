function rho_search()
%% general
i = 5;
x0 = [1, 800];
filename = 'inputHexane.mat';

dat = readtable('PRSV_data.xlsx');
Na = 6.02214076e23; % 2019 SI units precisely
k  = 1.380649e-23;  % 2019 SI units precisely
R  = Na*k;

name = char(dat{i,1});
Tc = dat{i,2};
pc = dat{i,3};
w  = dat{i,4};  
kappa_1 = dat{i,5};
kappa_2 = dat{i,6};
kappa_3 = dat{i,7};
TL = dat{i,8};

T  = [TL+0.15:3:Tc, Tc];

%% n_HEXANE
% x0 = [8, 7590];

%% n_HEPTANE
% x0 = [2, 6000];

%% n_OCTANE
% x0 = [2, 6000];

%% n_DECANE
% x0 = [2, 6000];

%% CYCLOHEXANE
% x0 = [4, 9200];

%% BENZENE
% x0 = [1, 11215];
%%

for i = 1:length(T)
    Tr = T(i)/Tc;
    kappa_0 = 0.378893 + 1.4897153*w - 0.17131848*w^2 + 0.0196544*w^3;
    if any(strcmp(name,{'Water','Ethanol','1-Decanol'}))
        if Tr >= 1
            kappa = kappa_0;
        else
            kappa   = kappa_0 + (kappa_1 + kappa_2*(kappa_3 - Tr)*(1 - sqrt(Tr)))*(1 + sqrt(Tr))*(0.7 - Tr);    %PRSV EoS
        end
    else
        if Tr >= 0.7
            kappa = kappa_0;
        else
            kappa   = kappa_0 + (kappa_1 + kappa_2*(kappa_3 - Tr)*(1 - sqrt(Tr)))*(1 + sqrt(Tr))*(0.7 - Tr);    %PRSV EoS
        end
    end
    %kappa = (0.37464+1.54226*w-0.26992*w^2);       % PR EoS
    alpha = (1+kappa*(1-sqrt(Tr)))^2;
    a = 0.457235*R^2*Tc^2/pc*alpha;
    b = 0.077796*R*Tc/pc;

    Ai = - 10.0e-16/(1.2326 + 1.3757*w);
    Bi = + 10.0e-16/(0.9051 + 1.5410*w);
    %di = Ai*(1 - Tr) + Bi;
    di = 0.305*Na^(-2/3);
    K = di*a*b^(2/3)/R/T(i);
    %rhoc = 1/vc;
    pars = [R,T(i),a,b];

    options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxIterations',1000,'MaxFunctionEvaluations',1000);
    %rho = fsolve(@(x) model(x,pars),x0,options)
    options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);
    %rho = fminsearch(@(x) model(x,pars),x0,options)
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
        'MaxFunctionEvaluations',10000000,'MaxIterations',10000000,'FunctionTolerance',1e-10,...
        'StepTolerance',1e-10);
    rho(i,:) = lsqnonlin(@(x) model(x,pars),x0,[0,0],[],options);

    mu1 = R*T(i)*log(abs(rho(i,1)/(1 - b*rho(i,1)))) + R*T(i)*b*rho(i,1)/(1 - b*rho(i,1)) ...
         + (a/2/sqrt(2)/b)*log(abs((1 + (1 - sqrt(2))*b*rho(i,1))/(1 + (1 + sqrt(2))*b*rho(i,1)))) ...
         - a*rho(i,1)/(1 + 2*b*rho(i,1) - b^2*rho(i,1)^2);
    mu2 = R*T(i)*log(abs(rho(i,2)/(1 - b*rho(i,2)))) + R*T(i)*b*rho(i,2)/(1 - b*rho(i,2)) ...
     + (a/2/sqrt(2)/b)*log(abs((1 + (1 - sqrt(2))*b*rho(i,2))/(1 + (1 + sqrt(2))*b*rho(i,2)))) ...
     - a*rho(i,2)/(1 + 2*b*rho(i,2) - b^2*rho(i,2)^2);
    
    
    f1 = R*T(i)*rho(i,1)*(log(abs(rho(i,1))) - 1) - R*T(i)*rho(i,1)*log(abs(1-b*rho(i,1))) ...
         + a*rho(i,1)/2/sqrt(2)/b*log(abs((1 + (1 - sqrt(2))*b*rho(i,1))/(1 + (1 + sqrt(2))*b*rho(i,1))));
    f2 = R*T(i)*rho(i,2)*(log(abs(rho(i,2))) - 1) - R*T(i)*rho(i,2)*log(abs(1-b*rho(i,2))) ...
         + a*rho(i,2)/2/sqrt(2)/b*log(abs((1 + (1 - sqrt(2))*b*rho(i,2))/(1 + (1 + sqrt(2))*b*rho(i,2))));
    p1 = -f1 + rho(i,1)*mu1;
    p2 = -f2 + rho(i,2)*mu2;
    
    input(i,:) = [T(i), rho(i,1), rho(i,2), mu1, mu2, p1, p2, Tr, K];
end

%disp(input)
save(filename,'input')
plot(T,rho(:,2))
legend('rhoL vs temp')

figure

plot(rho(:,1),input(:,6),'x',rho(:,2),input(:,7),'x')
legend('rhoV','rhoL')


%RHO = linspace(rho(1),rho(2),1000);
%RHO = linspace(0,7000,1000);

% for i=1:1000
%     MU(i) = R*T*log(abs(RHO(i)/(1 - b*RHO(i)))) + R*T*b*RHO(i)/(1 - b*RHO(i)) ...
%     + (a/2/sqrt(2)/b)*log(abs((1 + (1 - sqrt(2))*b*RHO(i))/(1 + (1 + sqrt(2))*b*RHO(i)))) ...
%     - a*RHO(i)/(1 + 2*b*RHO(i) - b^2*RHO(i)^2);
%     FF(i) = R*T*RHO(i)*(log(abs(RHO(i))) - 1) - R*T*RHO(i)*log(abs(1-b*RHO(i))) ...
%      + a*RHO(i)/2/sqrt(2)/b*log(abs((1 + (1 - sqrt(2))*b*RHO(i))/(1 + (1 + sqrt(2))*b*RHO(i))));
% end

% plot(RHO,MU)
% figure
% plot(1./RHO,FF)

end

function res = model(rho,pars)

R = pars(1);
T = pars(2);
a = pars(3);
b = pars(4);

mu1 = R*T*log(abs(rho(1)/(1 - b*rho(1)))) + R*T*b*rho(1)/(1 - b*rho(1)) ...
     + (a/2/sqrt(2)/b)*log(abs((1 + (1 - sqrt(2))*b*rho(1))/(1 + (1 + sqrt(2))*b*rho(1)))) ...
     - a*rho(1)/(1 + 2*b*rho(1) - b^2*rho(1)^2);
mu2 = R*T*log(abs(rho(2)/(1 - b*rho(2)))) + R*T*b*rho(2)/(1 - b*rho(2)) ...
 + (a/2/sqrt(2)/b)*log(abs((1 + (1 - sqrt(2))*b*rho(2))/(1 + (1 + sqrt(2))*b*rho(2)))) ...
 - a*rho(2)/(1 + 2*b*rho(2) - b^2*rho(2)^2);

f1 = R*T*rho(1)*(log(abs(rho(1))) - 1) - R*T*rho(1)*log(abs(1-b*rho(1))) ...
     + a*rho(1)/2/sqrt(2)/b*log(abs((1 + (1 - sqrt(2))*b*rho(1))/(1 + (1 + sqrt(2))*b*rho(1))));
f2 = R*T*rho(2)*(log(abs(rho(2))) - 1) - R*T*rho(2)*log(abs(1-b*rho(2))) ...
    + a*rho(2)/2/sqrt(2)/b*log(abs((1 + (1 - sqrt(2))*b*rho(2))/(1 + (1 + sqrt(2))*b*rho(2))));
p1 = -f1 + rho(1)*mu1;
p2 = -f2 + rho(2)*mu2;

res = [(mu1 - mu2); (p1 - p2)];
end

function res = fugacity()
end