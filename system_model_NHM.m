function [ dX ] = system_model_NHM(t,X)

global N J m alfa Ko Kd 

% Acting forces
[F1,F2,ang]=forces_NHM(X);
FT=F1+F2;

% Magnitude of F1
F_nV=(sqrt(sum(abs(F1).^2,2)));

% desired theta
thr=mod(ang',2*pi);

% actual theta
th=mod(X(3:6:end),2*pi);

ang=th-thr;
td=[ang ang+2*pi ang-2*pi];
[~,I]=min(abs(td),[],2);

dX=zeros(6*N,1);

for i=1:N
    a(i)=td(i,I(i));
    kl=0.3;
    kth=J(i)*kl*F_nV(i);
    kom=J(i)*(1+alfa)*sqrt(kl*F_nV(i)/alfa);
    dX(6*i-5)=X(6*i-2)*cos(X(6*i-3))-X(6*i-1)*sin(X(6*i-3));
    dX(6*i-4)=X(6*i-2)*sin(X(6*i-3))+X(6*i-1)*cos(X(6*i-3));
    dX(6*i-3)=X(6*i);
    dX(6*i-2)=1/m(i)*FT(i,:)*[cos(X(6*i-3)) sin(X(6*i-3))]';
    dX(6*i-1)=1/m(i)*(Ko*F2(i,:)*[-sin(X(6*i-3)) cos(X(6*i-3))]'-Kd*X(6*i-1));
    dX(6*i)  =1/J(i)*(-kth*a(i)-kom*X(6*i));
end

end

