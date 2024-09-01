clc
clear all
close all
global nameVectorVar; nameVectorVar='ke';
global vectorVar; vectorVar=[3,3.5,4,20];
global phitol; phitol=0.0001;
global phiend; phiend=0.6;
global zlim; zlim=10.0;
global ztol; ztol=0.01;
main();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = func_S(k,phi)
global vectorVar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global U1; U1 = 0;
global U2; U2=0;
global eta1; eta1=0.5;
global eta2; eta2=0.5;
global M; M=1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S=real((eta1*((M-U1)^2)*(1-(1-((2*phi)/(M-U1)^2))^0.5)) ...
+    (eta2*((M-U2)^2)*(1-(1-((2*phi)/(M-U2)^2))^0.5)) ...
+ (1-(1-((phi)/(vectorVar(k)-1.5)))^(1.5-vectorVar(k))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%% . NOTHING BEYOND THIS POINT SHOULD BE CHANGED







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main()
global nameVectorVar;
global vectorVar; 
for k=1:size(vectorVar,2)
[phiplot,Splot]=S_plot(k);
[z,phi,E] = phi_E(k);

hold on
figure(1)
plot(phiplot,Splot,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))))
xlabel('Phi')
ylabel('S')
title('S Plot')
legend show
% legend(Legend)


hold on
figure(2)
plot(z,phi,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))))
xlabel('Zeta')
ylabel('Phi')
title('Phi Plot')
legend show
% legend(Legend)


hold on
figure(3)
plot(z,E,'DisplayName',strcat(nameVectorVar,'= ', num2str(vectorVar(k))))
xlabel('Zeta')
ylabel('E')
title('E Plot')
legend show
% legend(Legend)


end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,phi,E] = phi_E(k)
global zlim;
global ztol;
global zero_point;
z=-zlim:ztol:zlim;
z0index=((size(z,2)-1)/2)+1;
phi=zeros(1,size(z,2));
phi(z0index)=zero_point;
for ii=z0index+1:size(phi,2)
phi(ii)=phi(ii-1)-(ztol)*  sqrt(-2*func_S(k,phi(ii-1)) );
end
for ii=z0index-1:-1:1
phi(ii)=phi(ii+1)-(ztol)*  sqrt(-2*func_S(k,phi(ii+1)) );
end
E=zeros(1,size(phi,2));
for ii=2:size(phi,2)-1
E(ii)=-(phi(ii+1)-phi(ii-1))/(2*ztol);
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