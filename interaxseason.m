function dMdt = interaxseason(t,M,m,bi,bj,Nr,Nd,b300,b301,b310,b311,f00,f01,f10,f11,alpha,beta,w,H,prod,remin,pfrag,dfrag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M(M<1E-9) = 0;
N = M./m(:);
m=m(:);
dM = zeros(size(N));
dMremin = zeros(size(dM));
dMfrag = zeros(size(dM));

% Aggregation


M(M<0) = 0;
M = M(:);
N = M(:)./m(:);
m=m(:);
dM = zeros(size(N));
dMremin = zeros(size(dM));
dMfrag = zeros(size(dM));

%% Aggregation


%N = N(:);

for k = 1:length(bi)
    
    ii = bi(k)+1;
    jj = bj(k)+1;
    mi = m(ii);
    mj = m(jj);
    mij = m(ii)+m(jj);
    d00 = b300(k) + 1;
    d01 = b301(k) + 1;
    d10 = b310(k) + 1;
    d11 = b311(k) + 1;
    
    
    if length(beta) ==1
        dN = alpha*beta*N(ii)*N(jj);
    else
        dN = alpha*beta(k)*N(ii)*N(jj);
    end
    

    if dN > 0
        dM(ii) = dM(ii)-dN*mi;
        dM(jj) = dM(jj)-dN*mj;
        dM(d00) = dM(d00) + f00(k)*dN*mij; 
        dM(d01) = dM(d01) + f01(k)*dN*mij; 
        dM(d10) = dM(d10) + f10(k)*dN*mij; 
        dM(d11) = dM(d11) + f11(k)*dN*mij;
    end
end

%% Degradation
dMremin = -remin.*M(:);

%% Fragmentation
fM = dfrag*pfrag.*M;
%dMfrag = - fM  + [fM(Nd+1:end); zeros(Nd,1)];
dMfrag = zeros(size(fM));

%% collecting

dMflux = -M.*w(:)./H;
dMdt = dM + prod(:)*(1-cos(2*pi*t/365))/2 + dMflux + dMremin + dMfrag;

end
