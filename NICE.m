function [Cov_NICE,Corr_NICE,L_NICE] = NICE(X,Y,fac)
Ne = size(X,2) ; 
[CorrXY,~] = corr(X',Y');
std_rho = (1-CorrXY.^2)/sqrt(Ne);
sig_rho = sqrt(sum(sum(std_rho.^2)));

expo2 = 2:2:8;
for kk = 1:length(expo2)
    L = abs(CorrXY).^expo2(kk);
    Corr_NICE = L.*CorrXY;
    if norm(Corr_NICE - CorrXY,'fro') > fac*sig_rho
        expo2 = expo2(kk);
        break
    end
end
expo1 = expo2-2;
rho_exp1 = CorrXY.^expo1;
rho_exp2 = CorrXY.^expo2;

al = 0.1:.1:1;
for kk=1:length(al)
    L = (1-al(kk))*rho_exp1+al(kk)*rho_exp2;
    Corr_NICE = L.*CorrXY;
    if kk>1 && norm(Corr_NICE - CorrXY,'fro') > fac*sig_rho
        Corr_NICE = PrevCorr;
        break
    elseif norm(Corr_NICE - CorrXY,'fro') > fac*sig_rho
        break
    end
    PrevCorr = Corr_NICE;
    L_NICE = L;
end
Vy = diag(std(Y,0,2));
Vx = diag(std(X,0,2));
Cov_NICE = Vx*Corr_NICE*Vy;
end