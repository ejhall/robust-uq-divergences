% load('SEsim_20170427.mat')
% `defaultpath` should point to the raw data in `uqii_data_2`
% i.e. defaultpath = '.../uqii_data_2/rawdata/';
% PRESCRIBED
%
% MESH (COARSE FEM)
nc = 2^6+1;
dxc = 1/(nc-1); % coarse mesh (for FEM computation)
xc = 0:dxc:1;
% SAMPLES
m = 100000; % max 500000
mSkp = 1000; % produce m/mSkp sample means (each of size mSkp)
mk = 0:mSkp:m;
%
%
nfig = 1;
for uu = 1:2
    switch uu
        case 1
            % NOMINAL MODEL 1a with parameters theta=(mu=0.8, v=4.0, s=0.005, nv=0.045)
            %   corresponds to model(1,5)
            outmod = 'mod1a';
            ne = 5;
            sindex = 1;
            nvindex = 5;
            ifs=1;
            sfun = @(j) (sindex+j); % s+epsilon corresponds to sindex+1
            nvfun = @(j) (nvindex); % nv index is fixed
            guessxp = [[10, 20, 30, 40, 50]; 
                [10, 20, 30, 40, 50];
                [100, 100, 100, 100, 100]];   % or [10, 25, 35, 45, 55]
            guessxm = [[5, 10, 15, 20, 25];
                [5, 10, 15, 20, 25];
                [10, 20, 30, 40, 50]];% NR guesses
            guessbp = [[10, 20, 30, 40, 50];
                [10, 20, 30, 40, 50];
                [10, 20, 30, 40, 50]]; %[1.5, 3.1, 5.7, 7.5, 8.7]];   % NR guesses
            guessbm = [[5, 15, 25, 35, 45];
                [5, 15, 25, 35, 45];
                [5, 15, 25, 35, 45]];  % NR guesses
            guesscp = [[10, 20, 30, 40, 50];
                [10, 20, 30, 40, 50];
                [10, 20, 30, 40, 50]];   % NR guesses
            guesscm = [[10, 20, 30, 40, 50];
                [10, 20, 30, 40, 50];
                [20, 30, 40, 50, 60]];  % NR guesses
            % x = linspace(1,50,2.5e2);
        case 2
            % NOMINAL MODEL 1b with parameters theta=(mu=0.8, v=4.0, s=0.005, nv=0.045)
            %   corresponds to model(1,5)
            outmod = 'mod1b';
            ne = 5;
            sindex = 1;
            nvindex = 5;
            ifs=0;
            sfun = @(j) (sindex);       % s index is fixed
            nvfun = @(j) (nvindex+j);   % nv+epsilon corresponds to nvindex+1
            % guessxp = [[20, 20, 20, 20, 20];
            %     [10, 10, 10, 10, 10];
            %     [50, 50, 50, 50, 50]];   % NR guesses 
            % guessxm = [[5, 5, 5, 5, 5];
            %     [10, 10, 10, 10, 10];
            %     [20, 20, 20, 20, 20]];% NR guesses
            % guessbp = [[25, 25, 25, 25, 25];
            %     [75, 75, 75, 75, 75];
            %     [100, 100, 100, 100, 100]];   % NR guesses
            % guessbm = [[25, 25, 25, 25, 25];
            %     [30, 30, 30, 30, 30];
            %     [8, 8, 8, 8, 8]];  % NR guesses
            % guesscp = [[30, 30, 30, 30, 30]; 
            %     [20, 20, 20, 20, 20];
            %     [10, 10, 10, 10, 10]];   % NR guesses
            % guesscm = [[5, 5, 5, 5, 5]; 
            %     [10, 10, 10, 10, 10];
            %     [10, 10, 10, 10, 10]];  % NR guesses
            guessxp = [[0.05, 0.10, 0.20, 0.25, 0.3];
                [0.01, 0.07, 0.11, 0.13, 0.15];
                [0.05, 0.10, 0.20, 0.25, 0.3]];   
            guessxm = [[0.05, 0.10, 0.20, 0.25, 0.3];
                [0.01, 0.06, 0.1, 0.11, 0.15];
                [0.05, 0.10, 0.20, 0.25, 0.3]];
            guessbp = [[0.05, 0.10, 0.20, 0.25, 0.3];
                [0.01, 0.07, 0.11, 0.13, 0.15];                       
                [0.05, 0.10, 0.20, 0.25, 0.3]];                       
            guessbm = [[0.05, 0.10, 0.20, 0.25, 0.3];
                [0.01, 0.07, 0.11, 0.13, 0.15];                        
                [0.05, 0.10, 0.20, 0.25, 0.3]];                       
            guesscp = [[0.05, 0.10, 0.20, 0.25, 0.3]; 
                [0.01, 0.07, 0.11, 0.13, 0.15];
                [0.05, 0.10, 0.20, 0.25, 0.3]];                       
            guesscm = [[0.05, 0.10, 0.2, 0.25, 0.3]; 
                [0.01, 0.07, 0.11, 0.13, 0.15];
                [0.05, 0.10, 0.20, 0.25, 0.3]];                      
            % x = linspace(0,75,2.5e2); % for choosing NR guess for tau

    end
    %        
    Pmod = model(sindex,nvindex);
    %SP = Pmod.Sigma;
    load(strcat(defaultpath,Pmod.dir,'sigs.mat'), 'SIGMA6'); %'-regexp', '[SIG].');
    SP = SIGMA6;
    T = cholcov(SP);
    SPinv = pinv(T)*pinv(T.');
    ncSP = size(SP,1);
    % PARAMETER VALUES TO COMPARE (FD)
    t = zeros(1,length(ne));
    for jj=1:ne
        if ifs==1
            t(jj) = model(sfun(jj),nvfun(jj)).s;
            t0 = Pmod.s;
        else 
            t(jj) = model(sfun(jj),nvfun(jj)).nv;
            t0 = Pmod.nv;
        end      
    end
    %
    epsilon = abs(t-t0);
    epsilonlab = cellstr(num2str(epsilon', 'epsilon=%0.4f'));
    %
    %
    % ENCODE different observables g{1}(u) for sol'n vec u;
    g1 = @(u) (u(:, end) > 1.2);
    g2 = @(u) (u(:,end) > 0.25 & u(:,end) < .75);
    g3 = @(u) min(u(:, end), ones(size(u(:,end)))*3);
    % g4 = @(u) u(:, end);
    % g5 = @(u) arrayfun(@(i) trapz(xc,u(i,:)), 1:size(u,1))';
    % g6 = @(u) 1./(1 + exp(-u(:,1)));
    g = {g1, g2, g3};
    glab = {'$\chi (u(1) > 1.2)$', '$\chi (.25< u(1) <.75)$', ...
        '$\min(u(1),3)$'};  %, 'u(1)', 'trapz(u)', 'logit(u(1))'};
    idx = [1,2,3];
    ng = length(idx);
    % ENCODE upper and lower bound for each observable 
    ub1 = 1;
    ub2 = 1;
    ub3 = 3;
    % ub4 = Inf;
    % ub5 = Inf;
    % ub6 = 1;
    ub = {ub1, ub2, ub3}; %, ub4, ub5, ub6};
    ublab = {'1', '1', '3'}; %, 'Inf', 'Inf', '1'};
    lb1 = 0;
    lb2 = 0;
    lb3 = 0;
    % lb4 = 0;
    % lb5 = 0;
    % lb6 = 0;
    lb = {lb1, lb2, lb3}; %, lb4, lb5, lb6};
    lblab = {'0', '0', '0'}; %, '0', '0', '0'};
    clear g1 g2 g3 g4 g5 ub1 ub2 ub3 ub4 ub5 ub6 lb1 lb2 lb3 lb4 lb5 lb6;
    %
    % DATA STRUCTURES
    %
    % RE (#epsilon by 1)-matrix
    RE = zeros(ne, 1);
    % SI (#g by #epsilon by #runs)-matrix sensitivty index
    % WE (#g by #epsilon by #runs)-matrix weak error
    SI = zeros(ng, ne, length(mk)-1);
    WE = zeros(ng, ne, length(mk)-1);
    %
    % XI (#g by #epsilon by #runs)-matrix
    XIp = zeros(ng, ne, length(mk)-1);
    XIm = zeros(ng, ne, length(mk)-1);
    % BI (#g by #epsilon by #runs)-matrix using Bennett inequality
    BIp = zeros(ng, ne, length(mk)-1);
    BIm = zeros(ng, ne, length(mk)-1);
    % CI (#g by #epsilon by #runs)-matrix using corollary to Bennett
    CIp = zeros(ng, ne, length(mk)-1);
    CIm = zeros(ng, ne, length(mk)-1);
    % OPTIMAL CONSTANTS (#g by #epsilon by #runs)-matrix of c* values
    cXIp = zeros(ng, ne, length(mk)-1);
    cXIm = zeros(ng, ne, length(mk)-1);
    cBIp = zeros(ng, ne, length(mk)-1);
    cBIm = zeros(ng, ne, length(mk)-1);
    cCIp = zeros(ng, ne, length(mk)-1);
    cCIm = zeros(ng, ne, length(mk)-1);
    % CALCULATIONS
    %
    % 0.   COMPUTE MEAN and VARIANCE for each observable
    Ep = zeros(ng, length(mk)-1);
    Vp = zeros(ng, length(mk)-1);
    % ELxp = zeros(ng,ne,length(mk)-1);
    % VLxp = zeros(ng,ne,length(mk)-1);
    % ELxm = zeros(ng,ne,length(mk)-1);
    % VLxm = zeros(ng,ne,length(mk)-1);
    for ii=1:ng
        for kk=1:(length(mk)-1)
        up = Pmod.u(mk(kk)+1:mk(kk+1),:);
        Ep(ii,kk) = mean(g{idx(ii)}(up));
        Vp(ii,kk) = var(g{idx(ii)}(up), 0); % normalizes by N-1 (for second moment about mean change 0 to 1)
        end;
    end;
    %
    % I.   Weak error for each observable, for each altJ, for mSkp samples
    for ii=1:ng
        for jj=1:ne
            for kk=1:(length(mk)-1)
                uq = model(sfun(jj), nvfun(jj)).u(mk(kk)+1:mk(kk+1),:);
                Eq = mean(g{idx(ii)}(uq));
                WE(ii,jj,kk) = Eq-Ep(ii,kk); % weak error
                SI(ii,jj,kk) = t0*(Eq-Ep(ii,kk))/epsilon(jj); % log deriv
            end;
        end;
    end;
    clear ii jj kk uq Eq;
    %
    % II.  UQII BOUNDS
    %
    % 1° COMPUTE RE matrix between Pmod and family of alternative models
    for jj=1:ne
        %SQ = model(sfun(j),nvfun(j)).Sigma;
        %SP = Pmod.Sigma;
        load(strcat(defaultpath,model(sfun(jj),nvfun(jj)).dir,'sigs.mat'), 'SIGMA6'); % '-regexp', '[SIG].');
        SQ = SIGMA6;
        %MQ = model(sfun(j),nvfun(j)).Mu;
        %MP = Pmod.Mu;
        MQ = ones(ncSP, 1)*model(sfun(jj),nvfun(jj)).mu;
        MP = ones(ncSP, 1)*Pmod.mu;
        RE(jj) = 0.5*(log(det(SP)/det(SQ)) + trace(SPinv*SQ) + (MQ-MP).'*SPinv*(MQ-MP) - ncSP);
    end
    clear i j SQ SP MQ MP T SPinv SIGMA6;
    % 
    % 2° COMPUTE UQ bounds
    % Gps = zeros(ng,ne,length(x));
    % Hps = zeros(ng,ne,length(x));
    % Lps = zeros(ng,ne,length(x));
    % Gms = zeros(ng,ne,length(x));
    % Hms = zeros(ng,ne,length(x));
    % Lms = zeros(ng,ne,length(x));
    for ii=1:ng
        for jj=1:ne
            for kk=1:(length(mk)-1)
                % COMPUTE XI BOUNDS::
                up = Pmod.u(mk(kk)+1:mk(kk+1),:);
                Y = g{idx(ii)}(up)-Ep(ii,kk);
                dM = mk(kk+1)-mk(kk);      
                % To find an 'optimal' c:
                % Let Gp(c) be d/dc [ Lambda(c;g)/c + RE(Q|P)/c ],
                % then zeros of Gp are min of [Lambda(c;g)/c + RE(Q|P)/c].
                % Similarly, for Gm but min of [Lambda(c,-g)/c ...]
                %
                % Gp = @(c) log(sum(exp(c*Y))/dM)/c + RE(j)/c;
                % Gm = @(c) log(sum(exp(c*-Y))/dM)/c + RE(j)/c;
                %          v-- derivative of outside      
                dGp = @(c) -log(sum(exp(c*Y))/dM) / c^2 +...
                    sum(Y.*exp(c*Y)) / sum(exp(c*Y)) / c  - RE(jj) / c^2;   
                    % ^-- derivative of inside             ^-- d/dc RE/c
                dGm = @(c) -log(sum(exp(c*-Y))/dM) / c^2 +...
                    sum(-Y.*exp(c*-Y)) / sum(exp(c*-Y)) / c - RE(jj) / c^2;
                % Gps(ii,jj,:) = arrayfun(dGp, x);
                % Gms(ii,jj,:) = arrayfun(dGm, x);
                %
                % Find zeros using Newton-Raphson method
                options = optimset('TolX',1e-10,'TolFun',1e-4); % set tol %TolX 1e-10, TolFun 1e-6
                [cp, resnorm, f, exitflag, output, jacob] = ...
                    newtonraphson(dGp, guessxp(idx(ii),jj), options);
                cXIp(ii,jj,kk) = cp;
                [cm, resnorm, f, exitflag, output, jacob] = ...
                    newtonraphson(dGm, guessxm(idx(ii),jj), options);
                cXIm(ii,jj,kk) = cm;
                clear resnorm f exitflag output jacob;  
                XIp(ii,jj,kk) = log(sum(exp(cp*Y)/dM))/cp + RE(jj)/cp;
                XIm(ii,jj,kk) = -log(sum(exp(cm*-Y)/dM))/cm - RE(jj)/cm;
                %                 ELxp(ii,jj,kk) = mean(exp(cp*Y));
                %                 VLxp(ii,jj,kk) = var(exp(cp*Y));
                %                 ELxm(ii,jj,kk) = mean(exp(cm*-Y));
                %                 VLxm(ii,jj,kk) = mean(exp(cm*-Y));
                %
                % COMPUTE BI (BENNETT) BOUNDS::
                options = optimset('TolX',1e-10,'TolFun',1e-4);
                Y = ub{idx(ii)}-Ep(ii,kk);
                Z = -lb{idx(ii)}+Ep(ii,kk);
                s2 = Vp(ii,kk);
                %
                hp = @(c) (Y^2*exp(-c*s2/Y) + s2*exp(c*Y)) / (Y^2+s2);         
                dHp = @(c) -log(hp(c))/c^2 + (s2*Y * (exp(c*Y) -...
                    exp(-c*s2/Y)))/c/hp(c) - RE(jj)/c^2;
                hm = @(c) (Z^2*exp(-c*s2/Z) + s2*exp(c*Z)) / (Z^2+s2);         
                dHm = @(c) -log(hm(c))/c^2 + (s2*Z * (exp(c*Z) -...
                    exp(-c*s2/Z)))/c/hm(c) - RE(jj)/c^2;
                % Hps(ii,jj,:) = arrayfun(dHp, x);
                % Hms(ii,jj,:) = arrayfun(dHm, x);
                %
                % Find zeros using Newton-Raphson method
                if uu==1
                    [cp, resnorm, f, exitflag, output, jacob] = ...
                        newtonraphson(dHp, guessbp(idx(ii),jj), options);
                else % uu==2
                    cp = guessbp(idx(ii),jj); % use c* approx by RE
                end
                cBIp(ii,jj,kk) = cp;               
                if uu==1
                    [cm, resnorm, f, exitflag, output, jacob] = ...
                        newtonraphson(dHm, guessbm(idx(ii),jj), options);
                else % uu==2
                    cm = guessbm(idx(ii),jj);
                end
                cBIm(ii,jj,kk) = cm;
                clear resnorm f exitflag output jacob;  
                BIp(ii,jj,kk) = log(hp(cp))/cp + RE(jj)/cp;
                BIm(ii,jj,kk) = -log(hm(cm))/cm - RE(jj)/cm;
                %
                % COMPUTE CI (BENNETT-a,b) BOUNDS::
                a = lb{idx(ii)} - Ep(ii,kk);
                b = ub{idx(ii)} - Ep(ii,kk);
                lp = @(c) (b*exp(c*a) - a*exp(c*b))/(b-a);
                dLp = @(c) -log(lp(c))/c^2 + a*b*(exp(c*a) -...
                    exp(c*b))/c/(b*exp(c*a)-a*exp(c*b)) - RE(jj)/c^2;
                lm = @(c) (-b*exp(c*-a) + a*exp(c*-b))/(-b+a);
                dLm = @(c) -log(lm(c))/c^2 + a*b*(exp(c*-a) -...
                    exp(c*-b))/c/(-b*exp(c*-a)+a*exp(c*-b)) - RE(jj)/c^2;
                % Lps(ii,jj,:) = arrayfun(dLp, x);
                % Lms(ii,jj,:) = arrayfun(dLm, x);
                %
                % Find zeros using Newton-Raphson method
                if uu==1
                    [cp, resnorm, f, exitflag, output, jacob] = ...
                        newtonraphson(dLp, guesscp(idx(ii),jj), options);
                else % uu==2
                    cp = guesscp(idx(ii),jj); % use c* approx by RE
                end
                cCIp(ii,jj,kk) = cp;
                if uu==1
                    [cm, resnorm, f, exitflag, output, jacob] = ...
                        newtonraphson(dLm, guesscm(idx(ii),jj), options);
                else % uu==2
                    cm = guesscm(idx(ii),jj); % use c* approx by RE
                end
                cCIm(ii,jj,kk) = cm;
                clear resnorm f exitflag output jacob;  
                CIp(ii,jj,kk) = log(lp(cp))/cp + RE(jj)/cp;
                CIm(ii,jj,kk) = -log(lm(cm))/cm - RE(jj)/cm;
            end
        end
    end
    clear ii jj kk;  
    
    % 
    % 3° POST PROCESS FOR R PLOTTING 
    %   (bound_id, g_id, epsilon_val, RE_val, c*_val, bound_val)
    names_bounds = cell(8,1);;
    r1 = ones(1500,1);
    r2 = vertcat(ones(500,1), 2*ones(500,1), 3*ones(500,1));
    r3a = ones(100,1) * epsilon;
    r3 = repmat(r3a(:),1,3);
    r3 = r3(:);
    r4a = ones(100,1) * RE';
    r4 = repmat(r4a(:),1,3);
    r4 = r4(:);
    %
    names_bounds{1} = 'WE';
    out1 = permute(WE, [3,2,1]);
    out1 = out1(:);
    mat1 = horzcat(r1,r2,r3,r4,r1.*0,out1);
    names_bounds{2} = 'SI';
    out2 = permute(SI, [3,2,1]);
    out2 = out2(:);
    mat2 = horzcat(r1.*2,r2,r3,r4,r1.*0,out2);
    names_bounds{3} = 'XIp';
    out3 = permute(XIp, [3,2,1]);
    out3 = out3(:);
    cut3 = permute(cXIp, [3,2,1]);
    cut3 = cut3(:);
    mat3 = horzcat(r1.*3,r2,r3,r4,cut3,out3);
    names_bounds{4} = 'XIm';
    out4 = permute(XIm, [3,2,1]);
    out4 = out4(:);
    cut4 = permute(cXIm, [3,2,1]);
    cut4 = cut4(:);
    mat4 = horzcat(r1.*4,r2,r3,r4,cut4,out4);
    names_bounds{5} = 'BIp';
    out5 = permute(BIp, [3,2,1]);
    out5 = out5(:);
    cut5 = permute(cBIp, [3,2,1]);
    cut5 = cut5(:);
    mat5 = horzcat(r1.*5,r2,r3,r4,cut5,out5);
    names_bounds{6} = 'BIm';
    out6 = permute(BIm, [3,2,1]);
    out6 = out6(:);
    cut6 = permute(cBIm, [3,2,1]);
    cut6 = cut6(:);
    mat6 = horzcat(r1.*6,r2,r3,r4,cut6,out6);
    names_bounds{7} = 'CIp';
    out7 = permute(CIp, [3,2,1]);
    out7 = out7(:);
    cut7 = permute(cCIp, [3,2,1]);
    cut7 = cut7(:);
    mat7 = horzcat(r1.*7,r2,r3,r4,cut7,out7);
    names_bounds{8} = 'CIm';
    out8 = permute(CIm, [3,2,1]);
    out8 = out8(:);
    cut8 = permute(cCIm, [3,2,1]);
    cut8 = cut8(:);
    mat8 = horzcat(r1.*8,r2,r3,r4,cut8,out8);

    bounds = vertcat(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8);
    
    save(strcat('20170728_', outmod, '.mat'), ...
        'bounds', 'names_bounds', 'glab', 'Ep', 'Vp', ...
        'ELxp', 'ELxm', 'VLxp', 'VLxm', '-v6');
    
    % % plot some XIp and weak errors
    % figure(nfig)
    % plot(reshape(cXIp(:,1,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('c*, epsilon_1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(reshape(XIp(:,1,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('Xi_p, epsilon_1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(reshape(WE(:,1,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('WE, epsilon_1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(reshape(cXIp(:,3,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('c*, epsilon_3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(reshape(XIp(:,3,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('Xi_p, epsilon_3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(reshape(WE(:,3,:),[3,100])');
    % legend('g1','g2','g3')
    % title(strcat('WE, epsilon_3, ',outmod))
    % nfig = nfig+1;
    
    % % plot the funcitons to be minimized, to aid in NR guess
    % figure(nfig)
    % plot(x, reshape(Gps(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gp1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Gps(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gp2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Gps(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gp3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Gms(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gm1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Gms(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gm2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Gms(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Gm3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hps(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hp1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hps(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hp2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hps(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hp3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hms(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hm1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hms(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hm2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Hms(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Hm3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lps(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lp1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lps(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lp2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lps(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lp3, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lms(1,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lm1, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lms(2,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lm2, ',outmod))
    % nfig = nfig+1;
    % 
    % figure(nfig)
    % plot(x, reshape(Lms(3,:,:),[ne,length(x)])')
    % legend({'eps_1','eps_2','eps_3','eps_4','eps_5'}, 'Location', 'southeast')
    % title(strcat('Lm3, ',outmod))
    % nfig = nfig+1;
end

% % plot Ep and Vp for observable
% figure(nfig)
% plot(Ep')
% legend('g1','g2','g3')
% title(strcat('Ep, ',outmod))
% nfig = nfig+1;
% 
% figure(nfig)
% plot(Vp')
% legend('g1','g2','g3')
% title(strcat('Vp, ',outmod))
% nfig = nfig+1;