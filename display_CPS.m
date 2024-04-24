function [] = display_CPS(x_hat, x_corr, D, a_hat, a_corr, p, q, fig_n, title_str)


%% Initialize
% grid
H = 10; %height of the grid (# cells)
L = 10; %length of the grid (# cells)
W = 100; %width of a square cell (cm)

room_grid = zeros(2,p);

for i=1:p
	room_grid(1,i) = floor(W/2)+ mod(i-1,L)*W; 
	room_grid(2,i) = floor(W/2)+ floor((i-1)/L)*W;
end



%% Construct vectors to plot
% find targets
targ = find(x_hat);
targ_corr = find(x_corr);

% find attacks
atk = find(a_hat);
atk_corr = find(a_corr);

% find sensors and attacks positions
sens = zeros(p,3);          % sens(i,:) = [#sensors,atk,atk_corr]
for i = 1:q
    [~,max_ind] = max(D(i,:));
    sens(max_ind,1) = sens(max_ind,1) + 1;

    % save if sensor under attack
    if ismember(i,atk)
        sens(max_ind,2) = 1;
    end
    if ismember(i,atk_corr)
        sens(max_ind,3) = 1;
    end
end
atk_pos = sens(:,2) > 0;          % to logical
atk_corr_pos = sens(:,3) > 0;     % to logical

% compute grids for sensor plotting and indexes of attacks
room_grid_sensor = zeros(2,q);
atk_index = false(q,1);
atk_corr_index = false(q,1);
j = 1;
for i = 1:p     % for every cell
    if sens(i) == 1             % if 1 sensor in cell
        room_grid_sensor(:,j) = room_grid(:,i);
        if atk_pos(i); atk_index(j) = true; end
        if atk_corr_pos(i); atk_corr_index(j) = true; end
        j = j + 1;
    elseif sens(i) == 2         % if 2 sensors in cell
        room_grid_sensor(:,j) = room_grid(:,i) + W/8;
        if atk_pos(i); atk_index(j) = true; end
        if atk_corr_pos(i); atk_corr_index(j) = true; end
        j = j + 1;
        room_grid_sensor(:,j) = room_grid(:,i) - W/8;
        if atk_pos(i); atk_index(j) = true; end
        if atk_corr_pos(i); atk_corr_index(j) = true; end
        j = j + 1;
    end
end



%% Plots
% plot targets
figure(fig_n)
plot(room_grid(1,x_corr), room_grid(2,x_corr),'s','MarkerSize',10, ...
    'MarkerEdgeColor',1/255*[40 208 220], 'MarkerFaceColor',1/255*[40 208 220], ...
    'DisplayName','Targets')
hold on
plot(room_grid(1,x_hat), room_grid(2,x_hat),'o','MarkerSize',14, ...
    'MarkerEdgeColor',1/255*[40 40 200],'LineWidth',1.2, ...
    'DisplayName','Estimated targets')

% plot sensors and mark under attack

% sensors
plot(room_grid_sensor(1,:), room_grid_sensor(2,:),'o','MarkerSize',6, ...
    'MarkerEdgeColor',1/255*[40 220 40],'MarkerFaceColor',1/255*[40 220 40], ...
    'DisplayName','Sensors')

% attacks
plot(room_grid_sensor(1,atk_corr_index), room_grid_sensor(2,atk_corr_index),'*','MarkerSize',10, ...
    'MarkerFaceColor',1/255*[220 40 40],'LineWidth',1.1, ...
    'DisplayName','Attacks')
plot(room_grid_sensor(1,atk_index), room_grid_sensor(2,atk_index),'o','MarkerSize',10, ...
    'MarkerEdgeColor',1/255*[220 40 200],'LineWidth',1.2,...
    'DisplayName','Estimated attacks')

% plot parameters
grid on
xticks(100:100:1000)
yticks(100:100:1000)
xlabel('(cm)')
ylabel('(cm)')
axis([0 1000 0 1000])
axis square



%% Prints
str = '';

% x_corr - targ_corr
if ~isempty(x_corr)
    str = sprintf('%s Target positions: \n ', str);
    for i = 1:length(targ_corr)
        str = sprintf('%s%d  ', str, targ_corr(i));
    end
end

% x_hat - targ
if ~isempty(x_hat)
    str = sprintf('%s\n\n Estimated positions: \n', str);
    for i = 1:length(targ)
        str = sprintf('%s%d  ', str, targ(i));
    end
end

% a_corr - atk_corr
if ~isempty(a_corr)
    str = sprintf('%s\n\n Estimated attacks: \n', str);
    for i = 1:length(atk_corr)
        str = sprintf('%s%d ', str, atk_corr(i));
    end
end

% a_hat - atk
if ~isempty(a_hat)
    str = sprintf('%s\n\n Sensors under attack: \n', str);
    for i = 1:length(atk)
        str = sprintf('%s%d ', str, atk(i));
    end
end

% display legend, side text and title
legend('Location', 'bestoutside')
text(1100,300,str);
title(title_str)

hold off


end