% `defaultpath` should point to the raw data in `uqii_data_2`
% i.e. defaultpath = '.../uqii_data_2/rawdata/';

%% 1. FORM FIM using DERIVATIVE of the COVARIANCE FUNCTION
%
% Compute the FIM based on the (analytical) derivative of the two-point
% covariance function and the corresponding covariance matrix SIGMAdtheta
% (previously computed in R using RandomFields).
%
% FIMU is for hyperparameters theta = (beta, sigma, ell, tau) such that 
%  mu = mu(beta) and SIMGA = SIGMA(theta) = SIGMA(sigma, ell, tau)
% FIMW is for hyperparameters theta = (beta, sigma^2, ell, tau^2).
%
% Relevant dataset for nominal model:
plab = {'mu', 'sigma', 'ell', 'tau'};

svec=[0.005, 0.045, 0.085]';
nvvec=[0.005, 0.045, 0.085]';
RU = zeros(length(svec)*length(nvvec), 4); % vec form of s by nv matrix
RW = zeros(length(svec)*length(nvvec), 4); % vec form of s by nv matrix
for i=1:3
    for j=1:3
        localpath = strcat('SEsim_mu=0.8_v=4_s=', num2str(svec(i)), '_nv=', num2str(nvvec(j)), '/');
        load(strcat(defaultpath, localpath, 'dsigs.mat'));
        load(strcat(defaultpath, localpath, 'params.mat'), 'dxf', 'mu');
        % First compute SIGMA**-1 using special Cholesky factorizaiton and pinv
        T = cholcov(SIGMA0);
        SIGinv = pinv(T)*pinv(T.');
        % Then form the FIM matrix
        b = ones(1/dxf+1,1)*mu;
        FIM_beta = b.'*SIGinv*b;
        U = zeros(3,3); % with respect to sigma and tau
        U(1,1) = 0.5*trace(SIGinv*SIGMAdsig*SIGinv*SIGMAdsig);
        U(1,2) = 0.5*trace(SIGinv*SIGMAdsig*SIGinv*SIGMAdell);
        U(1,3) = 0.5*trace(SIGinv*SIGMAdsig*SIGinv*SIGMAdtau);
        U(2,2) = 0.5*trace(SIGinv*SIGMAdell*SIGinv*SIGMAdell);
        U(2,3) = 0.5*trace(SIGinv*SIGMAdell*SIGinv*SIGMAdtau);
        U(3,3) = 0.5*trace(SIGinv*SIGMAdtau*SIGinv*SIGMAdtau);
        W = zeros(3,3); % with respect to sigma^2 and tau^2
        W(1,1) = 0.5*trace(SIGinv*SIGMAdv*SIGinv*SIGMAdv);
        W(1,2) = 0.5*trace(SIGinv*SIGMAdv*SIGinv*SIGMAdell);
        W(1,3) = 0.5*trace(SIGinv*SIGMAdv*SIGinv*SIGMAdnv);
        W(2,2) = 0.5*trace(SIGinv*SIGMAdell*SIGinv*SIGMAdell);
        W(2,3) = 0.5*trace(SIGinv*SIGMAdell*SIGinv*SIGMAdnv);
        W(3,3) = 0.5*trace(SIGinv*SIGMAdnv*SIGinv*SIGMAdnv);
        %
        FIMU_theta = U + triu(U,1)';
        FIMU = blkdiag(FIM_beta,FIMU_theta);
        FIMW_theta = W + triu(W,1)';
        FIMW = blkdiag(FIM_beta,FIMW_theta);
        %
        RU(i + (j-1)*length(svec), :) = sqrt(diag(FIMU));
        RW(i + (j-1)*length(svec), :) = sqrt(diag(FIMW));
        SU(i + (j-1)*length(svec), :) = sqrt(diag(FIMU))/sum(sqrt(diag(FIMU)));
        SW(i + (j-1)*length(svec), :) = sqrt(diag(FIMW))/sum(sqrt(diag(FIMW)));
        fprintf('position %1i corresponds to s = %1.4f and nv = %1.4f\n', i + (j-1)*length(svec), svec(i), nvvec(j))
    end;
end;

save('FIM_screening.mat', 'RU', 'RW', 'SU', 'SW', '-v6');

% %% 2. COMPUTE MEAN and VARIANCE for OBSERVABLE
% %
% % Compute the mean F and variance V of an observable g(u) for u previously
% % computed. 
% %
% % Relevant dataset:
% dataset = 'SEsim_PmodelSA_20170426.mat';
% load(strcat(defaultpathSet, dataset));
% % Prescribed:
% % Nnominal model with hyperparemters theta=(s=0.055, v=4.0, mu=0.8, nv=0.055)
% Pmod = model(1,1,1);
% % Mesh (coarse FEM)
% nc = 2^6+1;
% dxc = 1/(nc-1); % coarse mesh (for FEM computation)
% xc = 0:dxc:1;
% % Samples
% M = 100000; % total number of samples (number of samples to load)
% Mskp = M; % number of samples per run
% Mk = 0:Mskp:M;
% % Encode five different observables g{1}(u) for sol'n vec u 
% g1 = @(u) (u(:, end) > 3.1);
% g2 = @(u) (u(:,(nc-1)/2) > 1 & u(:,(nc-1)/2) > 2);
% g3 = @(u) min(u(:, end), ones(size(u(:,end)))*50);
% g4 = @(u) u(:, end);
% g5 = @(u) arrayfun(@(i) trapz(xc,u(i,:)), 1:size(u,1))';
% g = {g1, g2, g3, g4, g5};
% glab = {'\chi (u>3.1)', '\chi (1<u<2)', 'min(u(1),50)', 'u(1)', 'trapz(u)'};
% idx = [1];
% ng = length(idx);
% %
% F = zeros(ng, length(Mk)-1);
% V = zeros(ng, length(Mk)-1);
% for i=1:ng
%     for k=1:(length(Mk)-1)
%     up = Pmod.u(Mk(k)+1:Mk(k+1),:);
%     F(i,k) = mean(g{idx(i)}(up));
%     V(i,k) = var(g{idx(i)}(up), 0); % normalizes by N-1 (for second moment about mean change 0 to 1)
%     end;
% end;
% 
% clearvars -except FIMU FIMW g glab idx ng F V plab