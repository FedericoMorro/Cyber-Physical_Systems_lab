%% April 11, 2024 - S. M. Fosson
%% Dynamic CPS: 3 targets moving in a room
%% The dyanmics is given by matrix A
clc
close all
clear all
load('e03/localization.mat')
n = size(A,1);
q = size(D,1);
T = n+q;
O = [D, eye(q)];
for i=1:T
    O = [O ; D*A^i eye(q)];
end

%% Check the rank of O: O is not full rank...
size(O)
rank(O)

pause

%% Graphical representation of the dynamic CPS

H = 10; %height of the grid (# cells)
L = 10; %length of the grid (# cells)
W = 100; %width of a square cell (cm)

room_grid = zeros(2,n);


for i=1:n
	room_grid(1,i) = floor(W/2)+ mod(i-1,L)*W; 
	room_grid(2,i) = floor(W/2)+ floor((i-1)/L)*W;
end

xtrue = zeros(n,1);
support = randperm(n);
support = support(1:3); % I consider 3 targets
xtrue(support) = 1;
target = find(xtrue)

for move = 1:100
    xtrue = A*xtrue;
    target = find(xtrue)
    
    plot(room_grid(1,target), room_grid(2,target),'s','MarkerSize',9, 'MarkerEdgeColor',1/255*[40 208 220],'MarkerFaceColor',1/255*[40 208 220])
    grid on
    legend( 'Targets','Location','eastoutside')
   
    xticks(100:100:1000)
    yticks(100:100:1000)
    xlabel('(cm)')
    ylabel('(cm)')
    axis([0 1000 0 1000])
    axis square
    str = sprintf(' Time = %d', move);
    text(1100,900,str);
    pause(1)
    hold off
end


