clear all;
% PRESCRIBED
fprintf('\n\nFIXING parameters and data structures\n')
defaultpath = '/Users/erichall/Research/hybrid_div_rpde/backup_data/experiment3/run20170602b/';

k1=1;   % starting index corresponds ot X in aX.mat
k2=500; % ending index corresponds to X in aX.mat
sampperfile = 1000; % number of samples per file
M = (k2+1-k1)*sampperfile; % number of samples (number of samples to load)
%
nc = 220/2 + 1;
dxc = 1/(nc-1); % coarse mesh (for FEM computation)
nf = 220 + 1;
dxf = 1/(nf-1);
skp = dxc/dxf;
%
% experiment is to compare different geostatistical models where the
% covariance structure is fit from 'lacking' data

%
%% DATA STRUCTURE model
f1 = 'dir';                 % string:   path to sample and parameter files
v1 = '';
f2 = 'abar';                % vector:   driving path abar (element-wise)
v2 = zeros(M,nc);
f3 = 'Mu';                  % vector:   means (on coarse mesh)
v3 = zeros(nf,1);
f4 = 'Sigma';               % matrix:   covariance matrix (on coarse mesh)
v4 = zeros(nf,nf);
f5 = 'u';                   % vector:   sol'n on coarse mesh
v5 = zeros(M,nc);
f6 = 'iy';                  % vector:   indices sampled
v6 = [];
%
model(1:10) = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6);
clear f1 v1 f2 v2 f3 v3 f4 v4 f5 v5 f6 v6;
%
%% LOAD RAW DATA for each model into structure, compute abar
fprintf('\n\nLOAD RAW DATA for each model into structure, compute abar\n')
for i=1:10
    fprintf('Percent of true = %i\n', i*10)
        model(i).dir = strcat('EmpGauss_', num2str(i),'/');
        %
        aa = zeros(M,nf);
        for k=k1:k2
            load(strcat(defaultpath, model(i).dir, 'a', num2str(k), '.mat'));
            aa((k-1)*sampperfile+1:k*sampperfile,:) = a;
        end 
        % for each Ke, compute abar(Ke) = 1/|Ke| *(int_Ke a dx) by trapezoidal rule
        c = cumsum(aa, 2); % sum along rows
        cs = (c(:,skp+1:skp:end)-[zeros(M,1),c(:,skp:skp:end-2)]);
        ls = aa(:,1:skp:end-2); 
        rs = aa(:,skp+1:skp:end);
        model(i).abar = (cs - 0.5*ls - 0.5*rs)*dxf/dxc; %a1bar = (cs - 0.5*ls - 0.5*rs)*dxf/dx;
        % load parameters
        S = load(strcat(defaultpath, 'params2.mat'), ...
            strcat('SIGMA', num2str(i)), strcat('mu', num2str(i)));
        names = fieldnames(S);
        [sig, mu] = names{:};
        model(i).Mu = S.(mu);
        model(i).Sigma = S.(sig);
        Y = load(strcat(defaultpath, 'yis.mat'), strcat('iy',num2str(i)));
        names = fieldnames(Y);
        [idx] = names{:};
        model(i).iy = Y.(idx);
        clear aa c cs ls rs S names sig mu Y idx;

end;
clear defaultpath params i k k1 k2;

%% COMPUTE solution to rPDE
fprintf('\n\nCOMPUTE solution to rPDE\n')
for i=1:10
    fprintf('FEM calculations for  percent of true = %i ...\n', i*10)
    iStart = tic;
    u = zeros(M,nc);
    abar = model(i).abar;
    parfor m=1:M
        u(m,:) = FEM1D([0, 2200], nc, abar(m,:)', ones(nc,1), 0, 0.5)'; 
    end;
    model(i).u = u;
    fprintf('                   (loop took %1.4f s)\n', toc(iStart))
end;
clear i m abar u;

%% SAVE model
fprintf('\n\nSAVE model\n')
save -v7.3 /Users/erichall/Research/hybrid_div_rpde/data_sets/experiment3/EmpGauss_20170603_L.mat;
