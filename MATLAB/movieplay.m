% movie speed, 1=30 fps
step=10;


color=[];
cir=0:0.01:2*pi;
if N==n_groups(1) % if there are no groups random colors
    color=[rand(N,1) rand(N,1) rand(N,1)];
else %if there are groups a color for each group
    for i=1:length(n_groups)
        color=[color; ones(n_groups(i),1)*rand ones(n_groups(i),1)*rand ones(n_groups(i),1)*rand];
    end
end

%%% MOVIE
scrsz = get(0,'ScreenSize');
vid=VideoWriter('simulation');
open(vid);
h=figure('Color','w','Position',[1 1 scrsz(3) scrsz(4)]);
for tt=1:step:length(tspan)
    % Plot of the walls
    for i=1:num_walls
        plot(map_walls(2*i-1,:),map_walls(2*i,:),'k','LineWidth',2);
        hold on
    end
    
    % Plot of the pedestrians represented as circles
    for i=1:N
        plot(r(i)*cos(cir)+X(tt,6*i-5),r(i)*sin(cir)+X(tt,6*i-4),'Color',color(i,:),'LineWidth',2) % plot cerchi
        hold on
        plot(r(i)*cos(X(tt,6*i-3))+X(tt,6*i-5),r(i)*sin(X(tt,6*i-3))+X(tt,6*i-4),'ok')
        hold on
    end
    
    hold off
    axis equal
    axis off
    title('HSFM Simulation','Interpreter','Latex','FontSize',20)
    M=getframe(h);
    writeVideo(vid,M);
end
% close(vid)
% 
% % %%% Uncomment to rewatch the movie
% % h=figure('Color','w','Position',[1 1 scrsz(3) scrsz(4)]);
% % movie(h,M)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%% Trajectories
% figure('Color','w')
% hold on
% for i=1:num_walls
%     plot(map_walls(2*i-1,:),map_walls(2*i,:),'k');
% end
% axis equal 
% plot(X(:,1:4:end),X(:,2:4:end),'o');
% % plot(m_c(1,:),m_c(2,:),'*k')

% %%% MOVIE
% scrsz = get(0,'ScreenSize');
% vid=VideoWriter('2room_ASFM','MPEG-4');
% open(vid);
% h=figure('Color','w','Position',[1 1 scrsz(3) scrsz(4)]);
% for tt=1:step:length(tspan)
%     for i=1:num_walls
%         plot(map_walls(2*i-1,:),map_walls(2*i,:),'k','LineWidth',2);
%         hold on
%     end
%     for i=1:N
%         plot(r(i)*cos(cir)+X(tt,6*i-5),r(i)*sin(cir)+X(tt,6*i-4),'Color',color(i,:),'LineWidth',2) % plot cerchi
%         hold on
%         plot(r(i)*cos(X(tt,6*i-3))+X(tt,6*i-5),r(i)*sin(X(tt,6*i-3))+X(tt,6*i-4),'ok')
%         %quiver(X(tt,6*i-5),X(tt,6*i-4),X(tt,6*i-2)*cos(X(tt,6*i-3)),X(tt,6*i-2)*sin(X(tt,6*i-3)))
%         hold on
%     end
%     hold off
%     %axis([0 25 0 25])
%     axis equal
%     line([0 0],[-5 0],'Color','k','LineWidth',2)
%     line([7 7],[-5 7],'Color','k','LineWidth',2)
%     xlim([-6 13])
%     ylim([-5 15])
%     M=getframe(h);
%     writeVideo(vid,M);
% end
% close(vid)