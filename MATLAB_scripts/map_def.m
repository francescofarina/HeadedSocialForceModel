function [map_walls ] = map_def()
% MAP_DEF Create a map where the walls and the obstacles are represented as
% lines.
%
% OUTPUT:
% - map_walls: a matrix containing the line segments 
%
% The i-th segment, AB, is defined as a cell s{i}=[xa ya; xb yb]'.
%
% Insert all the segments

% Insert the semgents of the map 

% s{1} = [0     15;     15      15]';
% s{2} = [15    15;     15      3.5]';
% s{3} = [15    1.5;    15      0]';
% s{4} = [15    0;      0       0]';
% s{5} = [0     0;      0       15]';
% s{6} = [15    0;      25      0]';
% s{7} = [15    15;     25      15]';
% s{8} = [25    15;     25      13.5]';
% s{9} = [25    11.5;   25      0]';
% s{10}= [25    8;      30      8]';
% s{11}= [30    8;      30      15]';
% s{12}= [30    15;     28.5    15]';
% s{13}= [25    15;     27      15]';
% s{14}= [26    15;     26      23]';
% s{15}= [29    15;     29      23]';
% s{16}= [26    23;     26.5    23]';
% s{17}= [29    23;     28.5    23]';

% s{1}=[0 0; 0 14]';
% s{2}=[7 0; 7 10]';
% s{3}=[0 7.5; 2.5 7.5]';
% s{4}=[4.5 7.5; 7 7.5]';
% s{5}=[0 14; 11 14]';
% s{6}=[7 9; 14 9]';
% s{7}=[7 12; 7 14]';
% s{8}=[14 9; 14 14]';
% s{9}=[13 14; 14 14]';
% s{10}=[11 14; 11 19]';
% s{11}=[13 14; 13 19]';

% Corridor
s{1}=[0 0; 0 25]';
s{2}=[5 0; 5 25]';

n_s=length(s); % Number of segments

map_walls=[];
for i=1:n_s
    map_walls=[map_walls;s{i}];
end
        
end

