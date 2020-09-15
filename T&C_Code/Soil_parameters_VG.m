%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute soil parameters_VG    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Necessary correction for VG soil hydraulic parameters to get a finite
%%% value at the residual water content 
function[s_S,b]=Soil_parameters_VG(Phy1,Osat,Ohy,nVG,alpVG,GRAPH)
%%%INPUTS
%%% 1 Van-Genuchten, 1980  Par: alpha,[1/mm]  m, [-] n, [-] Osat [-], Ohy [-],  Ks [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_S = zeros(1,length(nVG)); b= zeros(1,length(nVG));
mVG= 1-1./nVG;
for i=1:length(nVG)
    %%%%%%%%
    Phy = Phy1*101.9368;
    %%%%%%%%
    P = @(Se) -(1./alpVG(i)).*((Se).^(-1./mVG(i))-1).^(1./nVG(i));
    Pd =@(Se) (1./Se.^(1/mVG(i)) - 1).^(1/nVG(i) - 1)./(alpVG(i)*mVG(i)*nVG(i)*Se.^(1/mVG(i) + 1));
    %%%%%
    G = @(Se,b) Phy*exp(Se.*-abs(b));
    Gd = @(Se,b) -Phy*abs(b)*exp(Se.*-abs(b));
    %%%%%%%%%%%%
    %X1: Se; X2: b
    F = @(X) [ ( (Pd(X(1)))-(Gd(X(1),X(2)))); ((P(X(1)))-(G(X(1),X(2)))) ];
    %%%%%%
    options = optimoptions('fsolve','Display','off','MaxIter',1e4,'MaxFunctionEvaluations',1e4);
    
    [X] = fsolve(F,[0.05,-0.1],options);
    
    s_S(i) = X(1); b(i)=X(2);
    
    if isreal(s_S(i))~=1 ||  isreal(b(i))~=1
        disp('FAILURE in the approximation of Van-Genuchten soil water retention curve')
        return
    end
end
%%%%%
i=GRAPH; %%% Figure only for one layer 
if GRAPH >= 1
    O=Ohy(i):0.001:Osat(i);
    Se = (O-Ohy(i))./(Osat(i)-Ohy(i));
    
    
    P1 = (1./alpVG(i)).*((Se).^(-1./mVG(i))-1).^(1./nVG(i)); %%% [mm]
    P2 = -Phy*exp(Se.*-abs(b(i)));
    P3 = (Se<s_S(i)).*P2 +(Se>=s_S(i)).*P1;
    
    col='-g';
    figure(2)
    subplot(1,1,1) ;
    set(gca,'FontSize',9);
    semilogy(O,P1*9.8/(10^6),col,'Linewidth',2.5)
    grid on ; hold on ;
    semilogy(O,P3*9.8/(10^6),'.r','Linewidth',2.5)
    %semilogy(O,h,'--g','Linewidth',2.5)
    xlabel('\theta [-]') ; ylabel('Water Potential \Psi [MPa]')%,'interpreter','latex');
    set(gca,'YDir','reverse')
end
return 





