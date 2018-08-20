clear all;
% PRESCRIBED
fprintf('\n\nFIXING parameters and data structures\n')
%defaultpath = '/Users/erichall/Documents/Backup_data/uqii_data_2/'; 
params = 'params.mat';
k1=1;   % starting index corresponds ot X in aX.mat
k2=100; % ending index corresponds to X in aX.mat
sampperfile = 1000; % number of samples per file
M = (k2+1-k1)*sampperfile; % number of samples (number of samples to load)
%
nc = 2^6+1;
dxc = 1/(nc-1); % coarse mesh (for FEM computation)
nf = 2^10+1;
dxf = 1/(nf-1);
skp = dxc/dxf;
%
% experiment is to vary correlation lengths, nugget variance, mean
%s = 0.005:0.01:0.1;
s= 0.005:0.01:0.095;
ms = length(s);
%nv = 0.005:0.01:0.1;
nv = 0.005:0.01:0.095;
mnv = length(nv);
% %mu = 0.8:0.1:1.3;
% mu = 0.8;

%
%% DATA STRUCTURE model
f1 = 'dir';                 % string:   path to sample and parameter files
v1 = '';
f2 = 's';                   % float:    s  = ell  "correlation length"
v2 = 0;
f3 = 'v';                   % float:    v  = sigma^2 "variance/sill"
v3 = 0; 
f4 = 'mu';                  % float:    mu "trend/mean"
v4 = 0;
f5 = 'nv';                  % float:    nv = tau^2 "nugget effect"
v5 = 0;
f6 = 'abar';
v6 = zeros(M,nc);
f7 = 'Mu';                  % vector:   means (on coarse mesh)
v7 = zeros(nc,1);
f8 = 'Sigma';               % matrix:   covariance matrix (on coarse mesh)
v8 = zeros(nc,nc);
f9 = 'u';
v9 = zeros(M,nc);
%
model(1:ms, 1:mnv) = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8,f9,v9);
clear f1 v1 f2 v2 f3 v3 f4 v4 f5 v5 f6 v6 f6b v6b f7 v7 f8 v8 f9 v9;
%
%% LOAD RAW DATA for each model into structure, compute abar
fprintf('\n\nLOAD RAW DATA for each model into structure, compute abar\n')
for i=1:ms
    fprintf('s = %1.f\n', i)
    for j=1:mnv
        fprintf('nv = %1.f\n', j)
        model(i,j).dir = strcat('SEsim_mu=0.8_v=4_s=', num2str(s(i)), '_nv=', num2str(nv(j)), '/');
        
        model(i,j).s = s(i);
        model(i,j).nv = nv(j);
        %
        aa = zeros(M,nf);
        for k=k1:k2
            load(strcat(defaultpath, model(i,j).dir, 'a', num2str(k), '.mat'));
            aa((k-1)*sampperfile+1:k*sampperfile,:) = a;
        end 
        % for each Ke, compute abar(Ke) = 1/|Ke| *(int_Ke a dx) by trapezoidal rule
        c = cumsum(aa, 2); % sum along rows
        cs = (c(:,skp+1:skp:end)-[zeros(M,1),c(:,skp:skp:end-2)]);
        ls = aa(:,1:skp:end-2); 
        rs = aa(:,skp+1:skp:end);
        model(i,j).abar = (cs - 0.5*ls - 0.5*rs)*dxf/dxc; %a1bar = (cs - 0.5*ls - 0.5*rs)*dxf/dx;
        % load parameters
        load(strcat(defaultpath, model(i,j).dir, params), 'SIGMA0', 'v', 'mu');
        model(i,j).mu = mu;
        model(i,j).Mu = ones(nc, 1)*mu;
        model(i,j).v = v;
        model(i,j).Sigma = SIGMA0;
        clear aa c cs ls rs;
    end;
end;
clear defaultpath params i j k k1 k2;

%% COMPUTE solution to rPDE
fprintf('\n\nCOMPUTE solution to rPDE\n')
for i=1:ms
    fprintf('FEM calculations for  s = %1.f ...\n', i)
    iStart = tic;
    for j=1:mnv
        fprintf('             ... for nv = %1.f ...', j)
        jStart = tic;
        u = zeros(M,nc);
        abar = model(i,j).abar;
        parfor m=1:M
        u(m,:) = FEM1D_abar(nc, abar(m,:)', ones(nc,1), 0, 1)'; 
        end;
        model(i,j).u = u;
        fprintf(' (nv loop took %1.4f s)\n', toc(jStart))
    end;
    fprintf('                             (s loop took %1.4f s)\n', toc(iStart))
end;
clear i j m abar u;

%% SAVE model
fprintf('\n\nSAVE model\n')
save -v7.3 uq_rpde/data_sets/SEsim_20170427.mat model;