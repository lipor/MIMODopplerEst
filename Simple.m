%% Simple Beampattern Matching Method

clear all

nT = 10;
L = 180;
R = zeros(nT);
theta = linspace(-90,90,L)*pi/180;
tt = theta*180/pi;
ppd = L/180; % points per degree
% phi = [zeros(60*ppd,1);ones(60*ppd,1);zeros(60*ppd,1)];
% phi = [zeros(30*ppd,1);ones(30*ppd,1);zeros(60*ppd,1);ones(30*ppd,1);zeros(30*ppd,1)];
% phi = [zeros(45*ppd,1);ones(2*ppd,1);zeros(86*ppd,1);ones(2*ppd,1);zeros(45*ppd,1)];
% phi = [zeros(45*ppd,1);ones(10*ppd,1);zeros(30*ppd,1);ones(10*ppd,1)...
%     ;zeros(30*ppd,1);ones(10*ppd,1);zeros(45*ppd,1)];
phi = [zeros(55*ppd,1);ones(10*ppd,1);zeros(20*ppd,1);ones(10*ppd,1)...
    ;zeros(45*ppd,1);ones(10*ppd,1);zeros(30*ppd,1)];

for i = 1:length(theta)
   aT = exp(  1j * pi*(0:nT-1)*sin(theta(i)) ).';
   R = R + phi(i)*aT*aT';
end
R = R/(nT);

P = zeros(L,1);
for i = 1:L
    aT = exp(  1j * pi*(0:nT-1)*sin(theta(i)) ).';
    P(i) = aT'*R*aT;
end
alpha = sum(phi)/sum(P);

% R2 = R2*sqrt(alpha);
% for i = 1:L
%     aT = exp(  1j * pi*(0:nT-1)*sin(theta(i)) ).';
%     aR = exp(  1j * pi*(0:nR-1)*sin(theta(i)) ).';
%     yt = kron(aT,aR);
%     P(i) = yt'*R2'*R2*yt; 
% end

plot(tt,10*log10(P/max(P)))