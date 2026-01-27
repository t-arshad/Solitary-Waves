clc
clear all
close all
global nameVectorVar; nameVectorVar='ke';
global vectorVar; vectorVar=[0.01];
global phitol; phitol=0.01;
global phiend; phiend=10.0;
global xilim; xilim=20.0;
global xitol; xitol=0.01;
global fontandaxessize;fontandaxessize=50;
global linewidth;linewidth=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('OCTAVE_VERSION', 'builtin') ~= 0
set(0, 'DefaultAxesFontSize', fontandaxessize)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = func_S(k,phi)
global vectorVar;
vv=vectorVar(k);

%//////////////////// This is the start of the function////////////////////

Ma=0.4;
Kz=0.5;
Kx=0.866;
zeta=1;
n=2;

r1=(Ma/Kz)^2;
r2=1-r1;
r3=3.0/(2.0*n);
r4=1.0/(2.0*n);
fn=(beta((zeta-r3),r3))/(beta((zeta-r4),r4));
cn=gamma(zeta)/(2*(gamma(1+r4))*(gamma(zeta-r4)));
a1=2*zeta*fn*cn*beta((zeta+r4),(1-r4));
a2=-(fn^2)*cn*beta((zeta+r3),(1-r3));
T1=-(0.5*a1)*(r2);
T2=(1.0/3.0)*((r1*a1^2)-(a2*r2));
S=(1.0/(Kx*Kx))*(T1*phi^2+T2*phi^3);

%//////////////////// This is the end of the function//////////////////////


%%%%%%%%%%%%%%%%%%%%%% Old - Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% . NOTHING BEYOND THIS POINT SHOULD BE CHANGED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main()
global linewidth
global nameVectorVar;
global vectorVar;
for k=1:size(vectorVar,2)
[phiplot,Splot]=S_plot(k);
[z,phi,E] = phi_E(k);

hold on
figure(1)
plot(phiplot,Splot,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
xlabel('\phi')
ylabel('V(\phi)')
title('V(\phi) Plot')
if(size(vectorVar,2)>1)
legend show
end

hold on
figure(2)
plot(z,phi,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
xlabel('\xi')
ylabel('\phi')
title('Phi Plot')
if(size(vectorVar,2)>1)
legend show
end

hold on
figure(3)
plot(z,E,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))),'linewidth',linewidth)
xlabel('\xi')
ylabel('E')
title('Energy Plot')
if(size(vectorVar,2)>1)
legend show
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
phi(ii)=phi(ii-1)-(xitol)*  sqrt(-2*func_S(k,phi(ii-1)) );
end
for ii=z0index-1:-1:1
phi(ii)=phi(ii+1)-(xitol)*  sqrt(-2*func_S(k,phi(ii+1)) );
end
E=zeros(1,size(phi,2));
for ii=2:size(phi,2)-1
E(ii)=-(phi(ii+1)-phi(ii-1))/(2*xitol);
end
E(1)=E(2);
E(end)=E(end-1);
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
S(ii) = func_S(k,phi(ii));
if S(ii)>0
revert_index=ii;
found_revert_index=1;
end
else
break;
end
end
phi=phi(1:size(S,2));
if revert_index==2
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
