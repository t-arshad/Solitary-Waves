clc
clear all
close all 


phitol=0.0001;
phiend=0.6;
zlim=10.0;
ztol=0.01;


U1=0;
U2=0;
eta1=0.5;
eta2=0.5;
M=1.1;
ke=[3,3.5,4,20];

for k=1:size(ke,2)
[phiplot,Splot,zero_point]=S_plot(U1,U2,eta1,eta2,M,ke(k),phitol,phiend);
[z,phi,E] = phi_E(zlim,ztol,zero_point,U1,U2,eta1,eta2,M,ke(k));

hold on
figure(1)
plot(phiplot,Splot,'DisplayName',strcat('ke=', num2str(ke(k))))
xlabel('Phi')
ylabel('S')
title('S Plot')
legend show
% legend(Legend)


hold on
figure(2)
plot(z,phi,'DisplayName',strcat('ke=', num2str(ke(k))))
xlabel('Zeta')
ylabel('Phi')
title('Phi Plot')
legend show
% legend(Legend)


hold on
figure(3)
plot(z,E,'DisplayName',strcat('ke=', num2str(ke(k))))
xlabel('Zeta')
ylabel('E')
title('E Plot')
legend show
% legend(Legend)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = func_S(U1,U2,eta1,eta2,M,ke,phi)
    S=real((eta1*((M-U1)^2)*(1-(1-((2*phi)/(M-U1)^2))^0.5)) ...
+    (eta2*((M-U2)^2)*(1-(1-((2*phi)/(M-U2)^2))^0.5)) ...
+ (1-(1-((phi)/(ke-1.5)))^(1.5-ke)));
end
% function S = func_S(U1,U2,eta1,eta2,M,ke, phi)
%     S=real((eta1*((M-U1)^2)*(1-(1-((2*phi)/(M-U1)^2))^0.5)) ...
% +    (eta2*((M-U2)^2)*(1-(1-((2*phi)/(M-U2)^2))^0.5)) ...
% + 1+3*(4*ke/(1+3*ke))-(1+3*(4*ke/(1+3*ke))*(1-phi)...
% +(4*ke/(1+3*ke))*(phi)^2)*exp(phi));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,phi,E] = phi_E(zlim,ztol,zero_point,U1,U2,eta1,eta2,M,ke)
z=-zlim:ztol:zlim;
z0index=((size(z,2)-1)/2)+1;
phi=zeros(1,size(z,2));
phi(z0index)=zero_point;
for ii=z0index+1:size(phi,2)
phi(ii)=phi(ii-1)-(ztol)*  sqrt(-2*func_S(U1,U2,eta1,eta2,M,ke,phi(ii-1)) );
end
for ii=z0index-1:-1:1
phi(ii)=phi(ii+1)-(ztol)*  sqrt(-2*func_S(U1,U2,eta1,eta2,M,ke,phi(ii+1)) );
end
E=zeros(1,size(phi,2));
for ii=2:size(phi,2)-1
E(ii)=-(phi(ii+1)-phi(ii-1))/(2*ztol);
end
E(1)=E(2);
E(end)=E(end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,S,zero_point] = S_plot(U1,U2,eta1,eta2,M,ke,phitol,phiend)
revert_index=-1;
found_revert_index=0;
phi=0.0:phitol:phiend;
for ii=1:size(phi,2)
if found_revert_index==0
S(ii) = func_S(U1,U2,eta1,eta2,M,ke,phi(ii));
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
disp(strcat('S IS POSITIVE FOR K=',num2str(ke))) 
elseif revert_index>-1
zero_point=phi(revert_index-1)+((0-S(revert_index-1))  /  ((S(revert_index)-S(revert_index-1))/(phi(revert_index)-phi(revert_index-1))));
phi(end)=zero_point;
S(end)=0;
else
zero_point=phi(end);
disp(strcat('ZERO POINT WAS NOT FOUND FOR K=',num2str(ke)))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%