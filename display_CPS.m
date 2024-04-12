function [] = display_CPS(x_hat, x, n, figure_n, str_title)

H = 10; %height of the grid (# cells)
L = 10; %length of the grid (# cells)
W = 100; %width of a square cell (cm)

room_grid = zeros(2,n);

for i=1:n
	room_grid(1,i) = floor(W/2)+ mod(i-1,L)*W; 
	room_grid(2,i) = floor(W/2)+ floor((i-1)/L)*W;
end

targ = find(x_hat);
targ_true = find(x);

figure(figure_n)
plot(room_grid(1,x_hat), room_grid(2,x_hat),'s','MarkerSize',9, 'MarkerEdgeColor',1/255*[40 208 220], 'MarkerFaceColor',1/255*[40 208 220])
hold on
plot(room_grid(1,x), room_grid(2,x),'s','MarkerSize',4, 'MarkerEdgeColor',1/255*[220 40 40], 'MarkerFaceColor',1/255*[220 40 40])
grid on

xticks(100:100:1000)
yticks(100:100:1000)
xlabel('(cm)')
ylabel('(cm)')
axis([0 1000 0 1000])
axis square

str = sprintf(' Estimated positions: \n');
for i=1:length(targ)
    str = sprintf('%s%d  ', str, targ(i));
end

if isempty(x)
    legend({'Estimated targets'}, 'Location', 'bestoutside')
else
    legend({'Estimated targets','Actual targets'}, 'Location', 'bestoutside')

    str = sprintf('%s\n\n Actual positions: \n ', str);
    for i=1:length(targ_true)
        str = sprintf('%s%d  ', str, targ_true(i));
    end
end

text(1100,500,str);
title(str_title)

hold off

end

