% load('EmpGauss_20170603_L.mat');
% PRESCRIBED
%
% SAMPLES
m = 500000; % max 500000
mSkp = 1000; % produce m/mSkp sample means (each of size mSkp)
mk = 0:mSkp:m;
% 
Pmod = model(5);
SP = Pmod.Sigma;
SPinv = pinv(SP);
MP = Pmod.Mu;
d = size(SP,1);
% NR guesses using sqrt(2*RE)/sqrt(Var_p(g))
guessxp = [[11.9204,    3.8359,    2.9665,    0.7417,   2.3792,    2.7231,    3.5818,    2.5299,    4.0320];
   [13.2122,    4.2516,    3.2880,    0.8220,    2.6371,   3.0182,    3.9699,    2.8040,    4.4689];
   [0,    0,    0,    0,    0,   0,    0,    0,    0];
   [26.3499,    8.4793,    6.5574,    1.6395,    5.2593,   6.0194,    7.9175,    5.5923,    8.9126]];
guessxm = [[11.9204,    3.8359,    2.9665,    0.7417,   2.3792,    2.7231,    3.5818,    2.5299,    4.0320];
   [13.2122,    4.2516,    3.2880,    0.8220,    2.6371,   3.0182,    3.9699,    2.8040,    4.4689];
   [0,    0,    0,    0,    0,   0,    0,    0,    0];
   [26.3499,    8.4793,    6.5574,    1.6395,    5.2593,   6.0194,    7.9175,    5.5923,    8.9126]];
% x = [linspace(0.1,1.9,2.5e1), linspace(2,11,1e1)];
%
% MODELS TO COMPARE TO
%
Qmods = [1, 2, 3, 4, 6, 7, 8, 9, 10];
nalt = length(Qmods);
Qps = 10*Qmods;
%
% ENCODE different observables g{1}(u) for sol'n vec u;
g1 = @(u) (u(:, end) > 6.2266e4);
g2 = @(u) (7.6256e4 > u(:,end) & u(:,end) > 4.8276e4);
g4 = @(u) u(:, end)/6.2266e4;
g = {g1, g2, g4};
glab = {'\chi (u(1) > 6.2266e4)', '\chi (7.6256e4 > u(1) > 4.8276e4)', ...
    'u(1)/6.2266e4'};
idx = [1,2,3];
ng = length(idx);
clear g1 g2 g4;
%
%% DATA STRUCTURES
%
% WE (#g by #altJ by #runs)-matrix
WE = zeros(ng, nalt, length(mk)-1);
% RE (#altJ by 1)-matrix
RE = zeros(nalt, 1);
% XI (#g by #altJ by #runs)-matrix
XIp = zeros(ng, nalt, length(mk)-1);
XIm = zeros(ng, nalt, length(mk)-1);
cXIp = zeros(ng, nalt, length(mk)-1);
cXIm = zeros(ng, nalt, length(mk)-1);

%% CALCULATIONS
%
% 0.   COMPUTE MEAN and VARIANCE for each observable
Ep = zeros(ng, length(mk)-1);
Vp = zeros(ng, length(mk)-1);
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
    for jj=1:nalt
        %W = zeros(length(mk)-1,1);
        for kk=1:(length(mk)-1)
            uq = model(Qmods(jj)).u(mk(kk)+1:mk(kk+1),:);
            Eq = mean(g{idx(ii)}(uq));
            WE(ii,jj,kk) = Eq-Ep(ii,kk);
        end;
    end;
end;
clear ii jj kk uq Eq W;
%
% II.  UQII BOUNDS
%
% 1° COMPUTE RE matrix between pMod and family of alternative models
% Gps = zeros(ng,nalt,length(x));
% Gms = zeros(ng,nalt,length(x));
for jj=1:nalt
    SQ = model(Qmods(jj)).Sigma;
    MQ = model(Qmods(jj)).Mu;
    %
    RE(jj) = 0.5 * (log(det(SP)) - log(det(SQ)) + ...
        trace(SPinv*SQ) + (MQ - MP)' * SPinv * (MQ - MP) - d);
end
clear ii jj mm SQ dSQ dSP detSP detSQ;
% 
% 2° COMPUTE XI bounds
for ii=1:ng
    for jj=1:nalt
        for kk=1:(length(mk)-1)
            up = Pmod.u(mk(kk)+1:mk(kk+1),:);
            Y = g{idx(ii)}(up)-Ep(ii,kk);
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
                sum(Y.*exp(c*Y)) / sum(exp(c*Y)) / c  - RE(jj) / c^2;   
                % ^-- derivative of inside             ^-- d/dc RE/c
            dGm = @(c) -log(sum(exp(c*-Y))/dM) / c^2 +...
                sum(-Y.*exp(c*-Y)) / sum(exp(c*-Y)) / c - RE(jj) / c^2;
            % Gps(ii,jj,:) = arrayfun(dGp, x);
            % Gms(ii,jj,:) = arrayfun(dGm, x);            
            %
            % Find zeros using Newton-Raphson method
            options = optimset('TolX',1e-6,'TolFun',1e-4); % set tol
            [cp, resnorm, f, exitflag, output, jacob] = ...
                newtonraphson(dGp, guessxp(idx(ii), jj), options);
            cXIp(ii,jj,kk) = cp;
            [cm, resnorm, f, exitflag, output, jacob] = ...
                newtonraphson(dGm, guessxm(idx(ii), jj), options);
            cXIm(ii,jj,kk) = cm;
            clear resnorm f exitflag output jacob;  
            XIp(ii,jj,kk) = log(sum(exp(cp*Y)/dM))/cp + RE(jj)/cp;
            XIm(ii,jj,kk) = -log(sum(exp(cm*-Y)/dM))/cm - RE(jj)/cm;
        end
    end
end
clear ii jj kk;
%

% PLOTS
%px=1;

% plot(Qmods, WE(:,:,1), 'b-x', Qmods, XIp(:,:,1), 'r-*', Qmods, XIm(:,:,1), 'r-*')
% hold on
% plot(Qmods, WE(:,:,1)+2*WE(:,:,2), 'b--', Qmods, WE(:,:,1)-2*WE(:,:,2), 'b--')

% % PLOT SENSITIVITY INDEX for deviation from nominal model
% f1 = figure(px);
% hold on
% W = WE(:,:,1);
% SE = reshape(permute(2*WE(:,:,2), [2 1]), nalt, 1, ng);
% [l,p] = boundedline(Qps, W, SE, '-g*', 'transparency', 0.1, 'alpha');
% Xp = XIp(:,:,1); % repmat(abs(epsilon),[3,1])
% XpE = reshape(permute(2*XIp(:,:,2), [2 1]), nalt, 1, ng); % reshape(repmat(abs(epsilon'),[3,1]), ne, 1, ng)
% Xm = XIm(:,:,1);
% XmE = reshape(permute(2*XIm(:,:,2), [2 1]), nalt, 1, ng);
% %
% [xpl,xpp] = boundedline(Qps, Xp, XpE, '-co', 'transparency', 0.1, 'alpha');
% [xml,xmp] = boundedline(Qps, Xm, XmE, '-ro', 'transparency', 0.1, 'alpha');
% % [bpl,bpp] = boundedline(abs(epsilon), Bp, BpE, '--gx', 'transparency', 0.1, 'alpha');
% % [bml,bmp] = boundedline(abs(epsilon), Bm, BmE, '--mx', 'transparency', 0.1, 'alpha');
% %
% outlinebounds(l,p);
% outlinebounds(xpl,xpp);
% outlinebounds(xml,xmp);
% %
% ax = gca;
% set(ax,'xlim',[Qps(1) Qps(end)]);
% %xlabel('Percent full data');
% %ylabel('Weak Error');
% hold off;
% px = px + 1;

save -v6 20170731_m5e5_mskp1e3_Pmod5_g124_outv6_R.mat W SE Xp Xm XpE XmE Qps;

% plot the funcitons to be minimized, to aid in NR guess
% nfig = 1;

% figure(nfig)
% plot(x, reshape(Gps(1,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gp1')
% nfig = nfig+1;
% 
% figure(nfig)
% plot(x, reshape(Gps(2,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gp2')
% nfig = nfig+1;
% 
% figure(nfig)
% plot(x, reshape(Gps(3,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gp3')
% nfig = nfig+1;
% 
% figure(nfig)
% plot(x, reshape(Gms(1,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gm1')
% nfig = nfig+1;
% 
% figure(nfig)
% plot(x, reshape(Gms(2,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gm2')
% nfig = nfig+1;
% 
% figure(nfig)
% plot(x, reshape(Gms(3,:,:),[nalt,length(x)])')
% legend({'Q_1','Q_2','Q_3','Q_4','Q_5','Q_6','Q_7','Q_8','Q_9','Q_{10}'}, 'Location', 'southeast')
% title('Gm3')
% nfig = nfig+1;