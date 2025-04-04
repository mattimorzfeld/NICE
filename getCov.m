function C = getCov(Nx, L, expType)
if expType == "singleScaleSpectral"
    C = spectralCov(Nx, L);
elseif expType == "singleScale"
    C = secov(Nx,L);
elseif expType == "exponential"
    C = ecov(Nx,L);
elseif expType == "multiScaleSpectral"
    C = L(2)*spectralCov(Nx, L(1)) + L(4)*spectralCov(Nx, L(3));
elseif expType == "multiScale"
    C = L(2)*secov(Nx,L(1)) + L(4)*secov(Nx,L(3));
elseif expType == "mcDSpectral"
    C = spectralCov(Nx, L);
    % derivative operator
    D = zeros(Nx);
    D = D+diag(ones(Nx-1,1),1)-diag(ones(Nx-1,1),-1);
    D(1, Nx) = -1;
    D(Nx, 1) = 1;
    D = 4*D;
    % form 2-variable matrix
    C = [C D*C; C*D' D*C*D'];
elseif expType == "mcD"
    C = secov(Nx,L);
    % derivative operator
    D = zeros(Nx);
    D = D+diag(ones(Nx-1,1),1)-diag(ones(Nx-1,1),-1);
    D(1, Nx) = -1;
    D(Nx, 1) = 1;
    D = 4*D;
    % form 2-variable matrix
    C = [C D*C; C*D' D*C*D'];
elseif expType == "NonStat"
    C = getCov_nonstationary(Nx);
elseif expType == "Wiggely"
    C = secovWig(Nx,L);
elseif expType == "Satellite"
    C = SatelliteCov(Nx,L);
elseif expType == "SpaceWeather"
   C = L(2)*secov(Nx,L(1)) + L(4)*secov(Nx,L(3));
else
    disp("Wrong expType = ", expType)
    C = zeros(Nx);
end
end

%% ----------------------------------------------------------
%% ----------------------------------------------------------
function C = getCov_nonstationary(n)

dx = 1;
a = .1;
b = 2;

C = zeros(n);
x = (1:n)*dx;
ell = a*x+b;

for ii=1:n
    for jj=ii:n
        dist = abs(ii-jj)*dx;
        lii = ell(ii);
        ljj = ell(jj);
        C(ii,jj) = (lii*ljj)^(1/4)/sqrt(mean([lii ljj])) * exp(-abs(dist/(mean([lii ljj])))^2);
    end
end

C = (C+C')-diag(diag(C));
end

function C = ecov(Nx,L)
C = zeros(Nx);
for ii=1:Nx
    for jj=ii:Nx
        d =min(abs(ii-jj), mod(-abs(ii-jj),Nx));
        C(ii,jj) = exp(-d/L);
    end
end
C = C+C'-diag(diag(C));
end
%% ----------------------------------------------------------
%% ----------------------------------------------------------

%% ----------------------------------------------------------
%% ----------------------------------------------------------
function C = secov(Nx,L)
C = zeros(Nx);
for ii=1:Nx
    for jj=ii:Nx
        d = min(abs(ii-jj), mod(-abs(ii-jj),Nx));
        C(ii,jj) = exp(-(d/L)^2);
    end
end
C = C+C'-diag(diag(C));
end
%% ----------------------------------------------------------
%% ----------------------------------------------------------

%% ----------------------------------------------------------
%% ----------------------------------------------------------
function C = secovWig(Nx,L)
C = zeros(Nx);
for ii=1:Nx
    for jj=ii:Nx
        d = min(abs(ii-jj), mod(-abs(ii-jj),Nx));
        C(ii,jj) = exp(-(d/L))*cos(d*Nx/pi);
    end
end
C = C+C'-diag(diag(C));
end
%% ----------------------------------------------------------
%% ----------------------------------------------------------

%% ----------------------------------------------------------
%% ----------------------------------------------------------
function C = spectralCov(Nx, L)
% create sinusoid matrix (basis functions)
% ... Nx must be even
% ... because we use sinusoids our correlation matrix will be periodic
E = zeros(Nx, Nx);
dz = 2*pi/Nx;
zj = (dz:dz:2*pi)';
E(:, 1) = ones(Nx, 1)/sqrt(Nx);
for i = 1:(Nx/2-1)
    tt = cos(i*zj);
    tn = sqrt(tt'*tt);
    E(:, 2*i) = tt/tn;
    tt = sin(i*zj);
    tn = sqrt(tt'*tt);
    E(:, 2*i+1) = tt/tn;
end
tt = cos(Nx*zj/2);
tn = sqrt(tt'*tt);
E(:, Nx) = tt/tn;
% create eigenvalues
b = 2*(pi*L/Nx)^2;
g = ones(Nx,1);
for i = 1:(Nx/2-1)
   g(2*i) = exp(-b*i^2);
   g(2*i+1) = g(2*i);
end
g(Nx) = exp(-b*(Nx/2)^2);
a = Nx/sum(g);
g = a*g;
G = diag(g);
C = E*G*E';
C = (C+C')/2;
end

function C = SatelliteCov(Nx,L)
C = zeros(Nx);
d1 = L(1);
d2 = L(2);
for ii=1:Nx
    for jj=ii:Nx
        C(ii,jj) = sqrt((ii*jj)/Nx^2)*exp(-0.5*((ii-jj)/d1)^2)...
                 + sqrt( (1-ii/Nx)*(1-jj/Nx) )*exp(-0.5*((ii-jj)/d2)^2);
    end
end
C = C+C'-diag(diag(C));
end
%% ----------------------------------------------------------
%% ----------------------------------------------------------
