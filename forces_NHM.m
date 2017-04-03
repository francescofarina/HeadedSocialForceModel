function [ F1,F2, ang ] = forces_NHM(X)

% FORCES determination of the resulting force on each individual
% Input X=[... x(i-1) y(i-1) dx(i-1) dy(i-1) x(i) y(i) dx(i-1) dy(i-1) ..]'

global n_groups
global N 
global r 
global m 
global v0  
global tau 
global A 
global B 
global Aw
global Bw

global k1 
global k2 
global e_act 
global e_ind 
global e_seq 
global e_n 
global num_walls 
global map_walls 

coef=0.5; % maximum distance from the waypoints we want to reach

% X decomposition
position=zeros(N,2);
vel=zeros(N,2);
for i=1:N
    position(i,:)=[X(6*i-5) X(6*i-4)];
    vel(i,:)=[X(6*i-2)*cos(X(6*i-3)) X(6*i-2)*sin(X(6*i-3))];
end


%Determination of the current waypoints
for i=1:n_groups(1)
    e(i,:)=(e_act{1}(i,:)'-position(i,:)')/norm(e_act{1}(i,:)'-position(i,:)');
    if norm(position(i,:)'-e_act{1}(i,:)')<=coef && e_ind{1}(i)<e_n{1}
        e_ind{1}(i)=e_ind{1}(i)+1;
        e_act{1}(i,:)=e_seq{1}(:,e_ind{1}(i));
        e(i,:)=(e_act{1}(i,:)'-position(i,:)')/norm(e_act{1}(i,:)'-position(i,:)');
    end
    if norm(position(i,:)'-e_act{1}(i,:)')<=coef && e_ind{1}(i)==e_n{1}
        e(i,:)=((1-exp(-5*norm(e_act{1}(i,:)'-position(i,:)')))/(1+exp(-5*norm(e_act{1}(i,:)'-position(i,:)'))))*((e_act{1}(i,:)'-position(i,:)')/norm(e_act{1}(i,:)'-position(i,:)'));
    end
end
for n=2:length(n_groups)
    for i=1:n_groups(n)
        e(sum(n_groups(1:n-1))+i,:)=(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)')/norm(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)');
        if norm(position(sum(n_groups(1:n-1))+i,:)'-e_act{n}(i,:)')<=coef && e_ind{n}(i)<e_n{n} %coef*r(sum(n_groups(1:n-1))+i)
            e_ind{n}(i)=e_ind{n}(i)+1;
            e_act{n}(i,:)=e_seq{n}(:,e_ind{n}(i));
            e(sum(n_groups(1:n-1))+i,:)=(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)')/norm(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)');
        end
        if norm(position(sum(n_groups(1:n-1))+i,:)'-e_act{n}(i,:)')<=coef && e_ind{n}(i)==e_n{n}
            e(sum(n_groups(1:n-1))+i,:)=((1-exp(-1*norm(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)')))/(1+exp(-1*norm(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)'))))*((e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)')/norm(e_act{n}(i,:)'-position(sum(n_groups(1:n-1))+i,:)'));
        end
    end
end

fi0=zeros(N,2); % velocity force
% Interindividual forces
fij1=zeros(N,2);% repulsive
fij2=zeros(N,2);% compression
fij3=zeros(N,2);% friction
% Obstacles
fiw1=zeros(N,2);% repulsive
fiw2=zeros(N,2);% compression
fiw3=zeros(N,2);% friction
% Groups
fgroup=zeros(N,2);
% e=([1 0]'*ones(1,N))';
for k=1:length(n_groups)
    for i=sum(n_groups(1:k))-n_groups(k)+1:sum(n_groups(1:k))
        fi0(i,:)=m(i)*(v0(i)*e(i,:)'-vel(i,:)')/tau;
        vect=(e(i,:)');
        ang(i)=atan2(vect(2),vect(1));
        for j=1:N
            if i~=j
                rij=r(i)+r(j);
                dij=norm(position(i,:)-position(j,:));
                nij=(position(i,:)'-position(j,:)')/dij;
                fij1(i,:)=fij1(i,:)'+A*exp((rij-dij)/B)*nij; 
                if dij<rij
                    fij2(i,:)=fij2(i,:)'+k1*(rij-dij)*nij;
                    tij=[-nij(2) nij(1)]';
                    dvij=(vel(j,:)'-vel(i,:)')'*tij;
                    fij3(i,:)=fij3(i,:)'+k2*(rij-dij)*dvij*tij;
                end 
            end
        end
        % Walls forces
        for w=1:num_walls
            xp=position(i,1);
            yp=position(i,2);
            rp=[xp yp]';
            ra=max([map_walls(2*w-1,1) map_walls(2*w,1)]',[map_walls(2*w-1,2) map_walls(2*w,2)]'); 
            rb=min([map_walls(2*w-1,1) map_walls(2*w,1)]',[map_walls(2*w-1,2) map_walls(2*w,2)]');
            xa=ra(1);
            ya=ra(2);
            xb=rb(1);
            yb=rb(2);
            % a point on AB can be parametrized as s(t)=ra+t(tb-ta), t in [0,1]
            % distance from s to p is phi(t)=||s(t)-p||
            % d(phi^2) gives the t which minimizes the distance from p to the
            % line in which AB lives. Since t in [0,1], t_star=min(max(0,t),1);
            % and the distance from p to AB is ||s(t_star)-p||
            t=((xp-xa)*(xb-xa)+(yp-ya)*(yb-ya))/(((xb-xa)^2+(yb-ya)^2));
            t_star=min(max(0,t),1);
            rh=ra+t_star*(rb-ra);
            diw=norm(rp-rh);
            niw=(rp-rh)/norm(rp-rh);
            tiw=[-niw(1) niw(2)]';
            fiw1(i,:)=fiw1(i,:)'+Aw*exp((r(i)-diw)/Bw)*niw;
            if diw<r(i)
                fiw2(i,:)=fiw2(i,:)'+k1*(r(i)-diw)*niw;
                fiw3(i,:)=fiw3(i,:)'-k2*(r(i)-diw)*(vel(i,:)*tiw)*tiw;
            end
        end
    end
end

% Force due to the desire to move as v0
F1=fi0;

% Other forces
F2=fij1+fij2+fij3+fiw1+fiw2+fiw3+fgroup;


end

