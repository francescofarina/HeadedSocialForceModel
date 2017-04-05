function [ dX ] = system_model_NHM(t,X)

global N J m alfa ko kd k1g k2g n_groups d_o d_f

% Positions and velocities
position=zeros(N,2);
vel=zeros(N,2);
for i=1:N
    position(i,:)=[X(6*i-5) X(6*i-4)];
    vel(i,:)=[X(6*i-2)*cos(X(6*i-3)) X(6*i-2)*sin(X(6*i-3))];
end

% Acting forces
[F0,Fe,ang]=forces_NHM(X);
FT=F0+Fe;

% Magnitude of F0
F_nV=(sqrt(sum(abs(F0).^2,2)));

% desired theta
thr=mod(ang',2*pi);

% actual theta
th=mod(X(3:6:end),2*pi);

% angle to rotate
ang=th-thr;
td=[ang ang+2*pi ang-2*pi];
[~,I]=min(abs(td),[],2);

dX=zeros(6*N,1);

% center o of mass of each group
ci=cell(1,length(n_groups));
for k=1:length(n_groups)
    ci{k}=[0 0]';
    for i=sum(n_groups(1:k))-n_groups(k)+1:sum(n_groups(1:k))
        ci{k}=ci{k}+position(i,:)';
    end
    ci{k}=ci{k}/n_groups(k);
end

for k=1:length(n_groups)
    for i=sum(n_groups(1:k))-n_groups(k)+1:sum(n_groups(1:k))
        a(i)=td(i,I(i));
        kl=0.3;
        kth=J(i)*kl*F_nV(i);
        kom=J(i)*(1+alfa)*sqrt(kl*F_nV(i)/alfa);

        p_i=(ci{k}-position(i,:)');

        dX(6*i-5)=X(6*i-2)*cos(X(6*i-3))-X(6*i-1)*sin(X(6*i-3));
        dX(6*i-4)=X(6*i-2)*sin(X(6*i-3))+X(6*i-1)*cos(X(6*i-3));
        dX(6*i-3)=X(6*i);

        % Here we substitute the step function in the definition of the group
        % cohesion forces with a sigmoid
        uf_group = k1g*(1+tanh(5*(abs(p_i'*[cos(X(6*i-3)) sin(X(6*i-3))]')-d_f)-3))*(p_i./norm(p_i))'*[cos(X(6*i-3)) sin(X(6*i-3))]';
        uo_group = k2g*(1+tanh(5*(abs(p_i'*[-sin(X(6*i-3)) cos(X(6*i-3))]')-d_o)-3))*(p_i./norm(p_i))'*[-sin(X(6*i-3)) cos(X(6*i-3))]';

        dX(6*i-2)=1/m(i)*(FT(i,:)*[cos(X(6*i-3)) sin(X(6*i-3))]'+uf_group);
        dX(6*i-1)=1/m(i)*(ko*Fe(i,:)*[-sin(X(6*i-3)) cos(X(6*i-3))]'-kd*X(6*i-1)+uo_group);
        dX(6*i)  =1/J(i)*(-kth*a(i)-kom*X(6*i));
    end
end

end

