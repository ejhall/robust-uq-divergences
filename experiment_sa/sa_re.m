% This file works off of the following data set:
% load('SEsim_20170427.mat')
% `defaultpath` should point to the raw data in `uqii_data_2`
% i.e. defaultpath = '.../uqii_data_2/rawdata/';
% PRESCRIBED
% experiment is to vary correlation lengths, nugget variance
s= 0.005:0.01:0.095;
ms = length(s);
nv= 0.005:0.01:0.095;
mnv = length(nv);

% FIX NOMINAL MODEL let model(P) be the nominal model
Pi = 1; % index of i in s
Pj = 5; % index of j in nv
P = Pi + (Pj-1)*ms; % i + (j-1)*length(s) for i in s, j in nv
Pmod = model(Pi,Pj);
load(strcat(defaultpath,model(Pi,Pj).dir,'sigs.mat'), 'SIGMA6');
%T = cholcov(Pmod.Sigma);
SP = SIGMA6;
T = cholcov(SP);
SPinv = pinv(T)*pinv(T.');
ncSP = size(SP,1);
%
% CALCULATE RELATIVE ENTROPY between nominal model and family of
% alternative models
RE = zeros(ms, mnv);
for i=1:ms
    for j=1:mnv
        %SQ = model(i, j).Sigma;
        %SP = Pmod.Sigma;
        load(strcat(defaultpath,model(i,j).dir,'sigs.mat'), 'SIGMA6'); %'-regexp', '[SIG].');
        SQ = SIGMA6;
        %MQ = model(i, j).Mu;
        %MP = Pmod.Mu;
        MQ = ones(ncSP, 1)*model(i,j).mu;
        MP = ones(ncSP, 1)*Pmod.mu;
        RE(i, j) = 0.5*(log(det(SP)/det(SQ)) + trace(SPinv*SQ) + (MQ-MP).'*SPinv*(MQ-MP) - ncSP);
    end;
end;
clear i j SQ SP MQ MP T SPinv;
RE
eps = s-Pmod.s;
epnv = nv-Pmod.nv;

% 
% PLOT RELATIVE ENTROPY LANDSCAPE
figure(4);
[X,Y] = meshgrid(epnv,eps);
surf(X,Y,RE);
hold all
shading interp;
dot = plot3(0,0,RE(P),'p','MarkerSize',14,'MarkerFaceColor','k','MarkerEdgeColor','w');
set(gca,'XTick', epnv); 
set(gca,'YTick', eps); 
xlabel('$\epsilon(\tau^2)$','Interpreter','LaTex');
ylabel('$\epsilon(\ell)$','Interpreter','LaTex');
zlabel('RE');
%title('RE landscape for (tau,ell) perturbations');
legend(dot,'Nominal model');
hold off;


figure(5)
hold all
contourf(epnv,eps,RE, 'LevelStep', 20.0);
dot = plot(0,0,'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
%colormap(summer);
xlabel('$\epsilon(\tau^2)$','Interpreter','LaTex');
ylabel('$\epsilon(\ell)$','Interpreter','LaTex');
set(gca,'XTick', epnv); 
set(gca,'YTick', eps); 
%title('RE landscape for $(\tau^2,\ell)$ perturbations','Interpreter','LaTex');
hcb = colorbar();
title(hcb,'RE')
legend(dot,'Nominal model');
hold off;
clear X Y;

% SAVE model
%save -v7.3 strcat('uqii_RE_', num2str(P),'.mat') model;

r3 = RE';
r3 = r3(:);
r2 = repmat(epnv',[10,1]);
r1 = repmat(eps,[10,1]);
r1 = r1(:);

% e(s), e(nv), RE
outmat = horzcat(r1,r2,r3);