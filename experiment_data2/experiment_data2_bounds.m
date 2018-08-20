% PRESCRIBED
%
% SAMPLES
m = 100000; % max 1e5
mSkp = 1000; % produce m/mSkp sample means (each of size mSkp)
mk = 0:mSkp:m;

for uu = 1:3
    switch uu
        case 1
            load('run20170608_a_v73.mat');
            outfile = '20170801a.mat';
            % NR guesses
            guess = [[2.2165,    2.0163,    0.9144];
                [0.4971,    0.4522,    0.2051];
                [0.5810,    0.5285,    0.2397];
                [0.5483,    0.4988,    0.2262]]';
            mu = 5.7548e+04;
            sdp = 6.9350e+04;
            sdm = 4.5746e+04;
            % x = [linspace(0.1,1.9,2.5e1), linspace(2,11,1e1)];
        case 2 
            load('run20170608_b_v73.mat');
            outfile = '20170801b.mat';
            % NR guesses
            guess = [[3.0499,    2.7745,    1.2582];
                [0.4972,    0.4523,    0.2051];
                [1.8801,    1.7103,    0.7756];
                [1.9278,    1.7537,    0.7953]]';
            mu = 6.8094e+04;
            sdp = 8.8179e+04;
            sdm = 4.8010e+04;
            % x = [linspace(0.1,1.9,2.5e1), linspace(2,11,1e1)];
        case 3
            load('run20170608_c_v73.mat');
            outfile = '20170801c.mat';
            % NR guesses
            guess = [[2.4973,    2.2717,    1.0302];
                [0.4972,    0.4523,    0.2051];
                [1.2798,    1.1642,    0.5280];
                [1.2571,    1.1435,    0.5186]]'; 
            mu = 6.8835e+04;
            sdp = 8.7010e+04;
            sdm = 5.0661e+04;
            % x = [linspace(0.1,1.9,2.5e1), linspace(2,11,1e1)];
    end
    %
    % NOMINAL MODEL 
    Pmod = model(1);
    % SP = Pmod.Sigma;
    % SPinv = pinv(SP);
    % MP = Pmod.Mu;
    % d = size(SP,1);
    %
    % MODELS TO COMPARE TO
    %
    allmodels = {'P', 'Qmax', 'Qmin', 'Qmean', 'Qmed'};
    altmodels = {'max', 'min', 'mean', 'med'};
    nalt = length(altmodels);
    %
    % ENCODE different observables g{1}(u) for sol'n vec u;
    g1 = @(u) (u(:, end) > mu);
    g2 = @(u) (sdp > u(:,end) & u(:,end) > sdm);
    g3 = @(u) u(:, end)/mu;
    g = {g1, g2, g3};
    glab = {'\chi (u(1) > mu)', '\chi (sdp > u(1) > sdm)', 'u(1)/mu'};
    idx = [1,2,3];
    ng = length(idx);
    clear g1 g2 g3;
    %
    % DATA STRUCTURES
    %
    % WE (#g by #altJ by 2)-matrix
    WE = zeros(ng, nalt, (length(mk)-1));
    % % RE (#altJ by 1)-matrix
    % RE = zeros(nalt, 1);
    % XI (#g by #altJ by #2)-matrix
    XIp = zeros(ng, nalt, (length(mk)-1));
    XIm = zeros(ng, nalt, (length(mk)-1));
    cXIp = zeros(ng, nalt, length(mk)-1);
    cXIm = zeros(ng, nalt, length(mk)-1);

    % CALCULATIONS
    %
    % % 0.   COMPUTE MEAN and VARIANCE for each observable
    % Ep = zeros(ng, length(mk)-1);
    % Vp = zeros(ng, length(mk)-1);
    % for ii=1:ng
    %     for kk=1:(length(mk)-1)
    %     up = Pmod.u(mk(kk)+1:mk(kk+1),:);
    %     Ep(ii,kk) = mean(g{idx(ii)}(up));
    %     Vp(ii,kk) = var(g{idx(ii)}(up), 0); % normalizes by N-1 (for second moment about mean change 0 to 1)
    %     end;
    % end;
    %
    % I.   Weak error for each observable, for each altmodel, for mSkp samples
    Ep = zeros(ng, nalt, length(mk)-1);
    Eq = zeros(ng, nalt, length(mk)-1);
    for ii=1:ng
        for jj=1:nalt
            for kk=1:(length(mk)-1)
                up = Pmod.u(mk(kk)+1:mk(kk+1),:);
                uq = model(jj+1).u(mk(kk)+1:mk(kk+1),:);
                Ep(ii,jj,kk) = mean(g{idx(ii)}(up)); %repetitive computation, but makes displaying UQ intervals easier
                Eq(ii,jj,kk) = mean(g{idx(ii)}(uq));
                WE(ii,jj,kk) = Eq(ii,jj,kk)-Ep(ii,jj,kk);
            end;
        end;
    end;
    clear ii jj kk uq up;
    %
    % II.  UQII BOUNDS
    %
    % 1° COMPUTE XI bounds
    % Gps = zeros(ng,nalt,length(x));
    % Gms = zeros(ng,nalt,length(x));
    for ii=1:ng
        for jj=1:nalt
            Xp = zeros(length(mk)-1,1);
            Xm = zeros(length(mk)-1,1);
            REQP = model(jj+1).RE
            for kk=1:(length(mk)-1)
                up = Pmod.u(mk(kk)+1:mk(kk+1),:);
                Y = g{idx(ii)}(up)-Ep(ii,jj,kk);
                dM = mk(kk+1)-mk(kk);      
                % To find an 'optimal' c:
                % Let Gp(c) be d/dc [ Lambda(c;g-F)/c + RE(Q|P)/c ],
                % then zeros of Gp are min of [Lambda(c;g-F)/c + RE(Q|P)/c].
                % Similarly, for Gm but min of [Lambda(c,-g+F)/c ...]
                %
                % Gp = @(c) log(sum(exp(c*Y))/dM)/c + RE(j)/c;
                % Gm = @(c) log(sum(exp(c*-Y))/dM)/c + RE(j)/c;
                %          v-- derivative of outside      
                dGp = @(c) -log(sum(exp(c*Y))/dM) / c^2 +...
                    sum(Y.*exp(c*Y)) / sum(exp(c*Y)) / c  - REQP / c^2;   
                    % ^-- derivative of inside             ^-- d/dc RE/c
                dGm = @(c) -log(sum(exp(c*-Y))/dM) / c^2 +...
                    sum(-Y.*exp(c*-Y)) / sum(exp(c*-Y)) / c - REQP / c^2;
                % Gps(ii,jj,:) = arrayfun(dGp, x);
                % Gms(ii,jj,:) = arrayfun(dGm, x);   
                %
                % Find zeros using Newton-Raphson method
                options = optimset('TolX',1e-10,'TolFun',1e-6); % set tol
                [cp, resnorm, f, exitflag, output, jacob] = ...
                    newtonraphson(dGp, guess(idx(ii),jj), options);

                [cm, resnorm, f, exitflag, output, jacob] = ...
                    newtonraphson(dGm, guess(idx(ii),jj), options);
                clear resnorm f exitflag output jacob;  
                cXIp(ii,jj,kk) = cp;
                cXIm(ii,jj,kk) = cm;
                XIp(ii,jj,kk) = log(sum(exp(cp*Y)/dM))/cp + REQP/cp;
                XIm(ii,jj,kk) = -log(sum(exp(cm*-Y)/dM))/cm - REQP/cm;
            end
        end
    end
    clear ii jj kk;

    save(outfile, 'WE', 'XIp', 'XIm', 'Eq', 'Ep');

end

% POST PROCESSIGN FOR R PLOTTING
load('20170731a.mat')
WEa = WE;
XIpa = XIp;
XIma = XIm;

load('20170731b.mat')
WEb = WE;
XIpb = XIp;
XImb = XIm;

load('20170731c.mat')
WEc = WE;
XIpc = XIp;
XImc = XIm;

r0a = ones(3600,1);
r0b = ones(3600,1).*2;
r0c = ones(3600,1).*3;

r1 = ones(1200,1);
r2 = vertcat( ones(400,1), 2*ones(400,1), 3*ones(400,1));
rr = ones(100,1) * (1:4);
r3 = repmat(rr(:),1,3);
r3 = r3(:);

outwa = permute(WEa, [3,2,1]);
outwa = outwa(:);
mat1a = horzcat(r1,r2,r3,outwa);

outpa = permute(XIpa,[3,2,1]);
outpa = outpa(:);
mat2a = horzcat(r1.*2, r2, r3, outpa);

outma = permute(XIma, [3,2,1]);
outma = outma(:);
mat3a = horzcat(r1.*3, r2, r3, outma);

outmata = horzcat(r0a, vertcat(mat1a, mat2a, mat3a));

%
outwb = permute(WEb, [3,2,1]);
outwb = outwb(:);
mat1b = horzcat(r1,r2,r3,outwb);

outpb = permute(XIpb,[3,2,1]);
outpb = outpb(:);
mat2b = horzcat(r1.*2, r2, r3, outpb);

outmb = permute(XImb, [3,2,1]);
outmb = outmb(:);
mat3b = horzcat(r1.*3, r2, r3, outmb);

outmatb = horzcat(r0b, vertcat(mat1b, mat2b, mat3b));

%
outwc = permute(WEc, [3,2,1]);
outwc = outwc(:);
mat1c = horzcat(r1,r2,r3,outwc);

outpc = permute(XIpc,[3,2,1]);
outpc = outpc(:);
mat2c = horzcat(r1.*2, r2, r3, outpc);

outmc = permute(XImc, [3,2,1]);
outmc = outmc(:);
mat3c = horzcat(r1.*3, r2, r3, outmc);

outmatc = horzcat(r0c, vertcat(mat1c, mat2c, mat3c));

%
%
outmat = vertcat(outmata, outmatb, outmatc);

save('20170731_idea2_outv6.mat', 'outmat', '-v6');