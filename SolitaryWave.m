clc
clear all
close all
global nameVectorVar; nameVectorVar='ke';
global vectorVar; vectorVar=[0.01];
global phitol; phitol=0.001;
global phiend; phiend=1e3;
global xilim; xilim=20.0;
global xitol; xitol=0.001;
global fontandaxessize;fontandaxessize=50;
global linewidth;linewidth=3;
global comparewithanalytical;comparewithanalytical=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ANY GLOBAL VARIABLES %%%%%%%%%%%%
global Ma;Ma=0.4;
global Kz;Kz=0.5;
global Kx;Kx=0.866;
global zeta;zeta=1;
global n;n=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('OCTAVE_VERSION', 'builtin') ~= 0
set(0, 'DefaultAxesFontSize', fontandaxessize)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V function starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = func_V(k,phi)
global vectorVar;
vv=vectorVar(k);
% Any parameters should be defined at the top as global variables %
%//////////////////// Current function starts///////////////////////////////
global Ma
global Kz
global Kx
global zeta
global n
r1=(Ma/Kz)^2;
r2=1-r1;
r3=3.0/(2.0*n);
r4=1.0/(2.0*n);
fn=(beta((zeta-r3),r3))/(beta((zeta-r4),r4));
cn=gamma(zeta)/(2*(gamma(1+r4))*(gamma(zeta-r4)));
a1=2*zeta*fn*cn*beta((zeta+r4),(1-r4));
a2=-(fn^2)*cn*zeta*beta((zeta+r3),(1-r3));
T1=-(0.5*a1)*(r2);
T2=(1.0/3.0)*((r1*a1^2)-(a2*r2));
S=(1.0/(Kx*Kx))*(T1*phi^2+T2*phi^3);
%//////////////////// Current function ends///////////////////////////////


%%%%%%%%%%%%%%%%%%%%%% Old - Functions starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U1=0;
% U2=0;
% eta1=0.5;
% eta2=0.5;
% M=1.1;
% S=real((eta1*((M-U1)^2)*(1-(1-((2*phi)/(M-U1)^2))^0.5)) ...
% +    (eta2*((M-U2)^2)*(1-(1-((2*phi)/(M-U2)^2))^0.5)) ...
% + (1-(1-((phi)/(vv-1.5)))^(1.5-vv)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U1=0;
%U2=0;
%eta1=0.5;
%eta2=0.5;
%M=1.4;
%S=real((eta1*((M-U1)^2)*(1-(1-((2*phi)/(M-U1)^2))^0.5)) ...
%+    (eta2*((M-U2)^2)*(1-(1-((2*phi)/(M-U2)^2))^0.5)) ...
%+ 1+3*(4*vv/(1+3*vv))-(1+3*(4*vv/(1+3*vv))*(1-phi)...
%+(4*vv/(1+3*vv))*(phi)^2)*exp(phi));

%%%%%%%%%%%%%%%%%%%%%% Old - Functions ends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V function ends %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analytic solution starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%//////////////////// Current analytic solution starts///////////////////////////////
function[V_analyt,E_analyt]=analytic_solution(xi_analyt)
global Ma
global Kz
global Kx
global zeta
global n
cn=gamma(zeta)./(2.*gamma(1+1./(2*n)).*gamma(zeta-1./(2*n)));
f=beta(zeta-3./(2*n),3./(2*n))./beta(zeta-1./(2*n),1./(2*n));
c1=2.*zeta.*f.*cn.*beta(zeta+1./(2*n),1-1./(2*n));
c2=-zeta.*f.^2.*cn.*beta(zeta+3./(2*n),1-3./(2*n));
f1=-(Ma^2.*c1)./Kz^2+c1;
f2=(Ma^2.*c2)./Kz^2+(Ma^2.*c1.^2)./Kz^2-c2;
A1=(2.*f2)./f1;
B1=Kx^2./f1;
p=3./A1;
w=sqrt(4.*B1);
dx=xi_analyt(2)-xi_analyt(1);
V_analyt=p.*(1./cosh(xi_analyt./w)).^2;
E_analyt=-gradient(V_analyt,dx);
end
%//////////////////// Current analytic solution ends///////////////////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analytic solution ends    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% . NOTHING BEYOND THIS POINT SHOULD BE CHANGED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main()
global linewidth
global nameVectorVar;
global vectorVar;
global comparewithanalytical;
kend=size(vectorVar,2);
if(comparewithanalytical)
kend=1;
end
for k=1:kend
[phiplot,Splot]=S_plot(k);
[z,phi,E] = phi_E(k);
if(comparewithanalytical)
[V_an,E_an]=analytic_solution(z);
end

figure(1)
hold on
plot(phiplot,Splot,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
xlabel('\phi')
ylabel('V(\phi)')
title('V(\phi) Plot')
if(size(vectorVar,2)>1&&!comparewithanalytical)
legend show
end

figure(2)
hold on
plot(z,phi,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
if (comparewithanalytical)
clf
hold on
plot(z,phi,'-o','DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
plot(z,V_an,'-x','DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
legend('Numerical','Analytical')
end
xlabel('\xi')
ylabel('\phi')
title('Phi Plot')
if(size(vectorVar,2)>1&&!comparewithanalytical)
legend show
end

figure(3)
hold on
plot(z,E,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
if (comparewithanalytical)
clf
hold on
plot(z,E,'-o','DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
plot(z,E_an,'-x','DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
legend('Numerical','Analytical')
end
xlabel('\xi')
ylabel('E')
title('Energy Plot')
if(size(vectorVar,2)>1&&!comparewithanalytical)
legend show
end

if (comparewithanalytical)
figure(4)
hold on
plot(z,abs(phi-V_an)./phi,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
xlabel('\xi')
ylabel('\phi Error')
title('Phi Error Plot')
disp(['Max error in φ: ' num2str(max(abs(phi-V_an)./phi)*100) ' %'])
if(size(vectorVar,2)>1&&!comparewithanalytical)
legend show
end

%figure(5)
%hold on
%plot(z,abs(E-E_an)./E,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
%xlabel('\xi')
%ylabel('E Error')
%title('Energy Error Plot')
%mask = E ~= 0;
%disp(['Max error in Energy: ' num2str(max(abs(E(mask)-E_an(mask))./E(mask))*100) ' %'])
%if(size(vectorVar,2)>1&&!comparewithanalytical)
%legend show
%end

end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,phi,E] = phi_E(k)
global xilim;
global xitol;
global zero_point;
z=-xilim:xitol:xilim;
z0index=((size(z,2)-1)/2)+1;
phi=zeros(1,size(z,2));
phi(z0index)=zero_point;
for ii=z0index+1:size(phi,2)
phi(ii)=phi(ii-1)-(xitol)*  sqrt(-2*func_V(k,phi(ii-1)) );
end
for ii=z0index-1:-1:1
phi(ii)=phi(ii+1)-(xitol)*  sqrt(-2*func_V(k,phi(ii+1)) );
end
E=zeros(1,size(phi,2));
for ii=2:size(phi,2)-1
E(ii)=-(phi(ii+1)-phi(ii-1))/(2*xitol);
end
E(1)   = 0;
E(end) = E(end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,S] = S_plot(k)
global phitol;
global phiend;
global zero_point;zero_point=-1;
global nameVectorVar;
global vectorVar;
revert_index=-1;
found_revert_index=0;
phi=0.0:phitol:phiend;
for ii=1:size(phi,2)
if found_revert_index==0
S(ii) = func_V(k,phi(ii));
if S(ii)>0
revert_index=ii;
found_revert_index=1;
end
else
break;
end
end
phi=phi(1:size(S,2));
if (revert_index==1)
zero_point=phi(end);
disp(strcat('S IS POSITIVE FOR ->',nameVectorVar,'=',num2str(vectorVar(k))))
elseif revert_index>-1
zero_point=phi(revert_index-1)+((0-S(revert_index-1))  /  ((S(revert_index)-S(revert_index-1))/(phi(revert_index)-phi(revert_index-1))));
phi(end)=zero_point;
S(end)=0;
else
zero_point=phi(end);
disp(strcat('ZERO POINT WAS NOT FOUND FOR ->',nameVectorVar,'=',num2str(vectorVar(k))))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main(); %Defined here because Octave is stupid
