clc; clearvars; close all; format short g; format compact;
tfinal = 30;
pars.D = 0.1; pars.Sf = 5; pars.mu = 0.5;
pars.K = 0.25;
% pars.K = 0.005;
pars.Yxs = 0.75; pars.Yps = 0.65;
D = pars.D; Sf = pars.Sf; mu = pars.mu;
K = pars.K; Yxs = pars.Yxs; Yps = pars.Yps;
cO = [5 0.02 0];
[t,c] = ode23s(@(t,c) chemofun (t,c,pars), [0 tfinal],cO);
nt = numel(t);
Sss = c(nt,1); Xss = c(nt,2); Pss = c(nt,3);
figure('units', 'normalized', 'outerposition', [0 0 1 1])
plot(t, c, 'LineWidth', 2); grid on;
xlabel('time, s', 'FontSize', 14); ylabel('concentration', 'FontSize', 14); legend('S', 'X', 'P');
clear t c
% Two steady states
Sss1 = D*K/(mu*Yxs-D); Xss1 = (Sf-Sss1)*Yxs;
Sss2 = Sf; Xss2 = 0;
state = 1
if state == 1
    Sss = Sss1; Xss = Xss1;
else
    Sss = Sss2; Xss = Xss2;
end
css = [Sss Xss]

A = [-D - mu*Xss*K/(K+Sss)^2, -mu*Sss/(K+Sss); mu*Xss*K*Yxs/(K+Sss)^2, -D+mu*Sss*Yxs/(K+Sss)];
[evec, eval] = eig(A)
evec1 = [evec(1,1); evec(2,1)]; evec2 = [evec(1,2); evec(2,2)];


% Phase plane for linearized equations
figure('units', 'normalized', 'outerposition',[0 0 1 1])
plot (Sss, Xss, 'ko', 'MarkerSize', 18, 'MarkerFaceColor', 'black');
hold on
% SO = 4; XO = 0.2;
% SOvec = 50; XOvec = X);
if state == 1
    SOvec = [0.95*Sss : 0.02*Sss: 1.05*Sss]; XOvec = [0.95*Xss: 0.02*Xss: 1.05*Xss];
else
    SOvec = [4.5 :0.1: 5.5]; XOvec = [0: 0.025: 0.1];
end
for i = 1:numel(SOvec)
    for j = 1:numel(XOvec)
        SO = SOvec(i);
        XO = XOvec(j);
        [t,c] = ode15s(@(t,c) chemofun2(t, c, pars, A, css), [0 tfinal], [SO XO]);
        % figure (2), plot(t,c)
        % pause
        Svec = c(:,1); Xvec = c(:,2);
        plot (Svec, Xvec, 'o-')
        xlabel('Substrate concentration', 'FontSize', 14); ylabel('Cell concentration', 'FontSize', 14);
        hold on;
        plot (SO ,XO, 'o', 'MarkerSize',7, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', [1 0.6 0.6]);
        hold on;
        % pause
    end
end

function f = chemofun (t, c, pars)
    S = c(1); X = c(2); P = c(3);
    D = pars.D; Sf = pars.Sf; mu = pars.mu; K = pars.K; Yxs = pars.Yxs; Yps = pars.Yps;
    rg = mu*(S/(K+S))*X;
    f(1,1) = D*(Sf-S)-rg;
    f(2,1) = -D*X + rg*Yxs;
    f(3,1) = -D*P + rg*Yps;
end
function f = chemofun2(t, c, pars, A, css)
    S = c(1); X = c(2);
    D = pars.D; Sf = pars.Sf;mu = pars.mu; K = pars.K; Yxs = pars. Yxs;
    Sss = css (1); Xss = css (2);
    f(1,1) = A(1,1)*(S-Sss) + A(1,2)*(X-Xss);
    f(2,1) = A(2,1)*(S-Sss) + A(2,2)*(X-Xss);
end