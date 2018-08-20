clear all;
% PRESCRIBED
fprintf('\n\nFIXING parameters and data structures\n')
defaultpath = '/Users/erichall/Research/hybrid_div_rpde/backup_data/experiment3/run20170608_c/';

k1=1;   % starting index corresponds ot X in aX.mat
k2=100; % ending index corresponds to X in aX.mat
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
allmodels = {'P', 'Qmax', 'Qmin', 'Qmean', 'Qmed'};
altmodels = {'max', 'min', 'mean', 'med'};
%model
%% DATA STRUCTURE model
f1 = 'dir';                 % string:   path to sample and parameter files
v1 = '';
f2 = 'abar';                % vector:   driving path abar (element-wise)
v2 = zeros(M,nc);
f3 = 'RE';                  % scalar:   means (on coarse mesh)
v3 = 0;
f4 = 'u';                   % vector:   sol'n on coarse mesh
v4 = zeros(M,nc);
f5 = 'iy';                  % vector:   indices sampled
v5 = [];
%
model(1:length(allmodels)) = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);
clear f1 v1 f2 v2 f3 v3 f4 v4 f5 v5;
%
%% LOAD RAW DATA for each model into structure, compute abar
fprintf('\n\nLOAD RAW DATA for each model into structure, compute abar\n')
for i=1:length(allmodels)
    fprintf('Loading data for model(%s)\n', allmodels{i})
    model(i).dir = strcat('Gauss_', allmodels{i},'/');
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
    clear aa c cs ls rs;
end
% load RE
load(strcat(defaultpath, 'RE.mat'));
model(1).RE = 0;
model(1).iy = yP;
itr = [imax, imin, imean, imed];
for i=1:4
    model(i+1).RE = RE(itr(i));
    model(i+1).iy = yQ(itr(i));
end
clear defaultpath itr i k k1 k2 RE imax imin imean imed yQ yP;

%% COMPUTE solution to rPDE
fprintf('\n\nCOMPUTE solution to rPDE\n')
for i=1:length(allmodels)
    fprintf('FEM calculations for  model(%s)...\n', allmodels{i})
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
save -v7.3 /Users/erichall/Research/hybrid_div_rpde/data_sets/experiment3/idea2/run20170608_c_v73.mat;
save -v6 /Users/erichall/Research/hybrid_div_rpde/data_sets/experiment3/idea2/run20170608_c_v6.mat;
