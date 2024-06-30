function [Adj, g] = network_config(type, N, silent)

Adj = zeros(N,N);
g = zeros(N,1);

if type == "linear"
    % linear configuration
    Adj(2,1) = 1;
    Adj(3,2) = 1;
    Adj(4,3) = 1;
    Adj(5,4) = 1;
    Adj(6,5) = 1;
    g(1) = 1;

elseif type == "tree"
    % binary tree configuration
    Adj(3,1) = 1;
    Adj(5,1) = 1;
    Adj(4,2) = 1;
    Adj(6,2) = 1;
    g(1) = 1;
    g(2) = 1;

else
    % selected configuration
    Adj(1,3) = 1;
    Adj(2,6) = 1;
    Adj(3,1) = 1;
    Adj(3,4) = 1;
    Adj(4,1) = 1;
    Adj(4,5) = 1;
    Adj(5,2) = 1;
    Adj(5,4) = 1;
    Adj(6,5) = 1;
    g(1) = 1;
    g(2) = 1;

end

aug_graph = [
    zeros(1,N+1)
    g       Adj
];
digr = digraph(aug_graph', {'0','1','2','3','4','5','6'});

if ~silent
    f = figure();
    f.Position([3 4]) = [525, 400];
    plot(digr, 'EdgeLabel', digr.Edges.Weight)
    title('Communication Network Topology')
end

end