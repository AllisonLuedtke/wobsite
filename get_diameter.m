%This function determines the maximum shortest path length for a network.
%It takes as input an adjacency matrix, Adj
%It returns the maximum shortest path length, max_shortest_path

%This uses Joseph Kirk's "dijkstra.m" path finding algorithm.

function [max_shortest_path] = get_diameter(Adj)

%Get the function variables
temp_A = Adj; %Ratchet adjacency matrix

%How many firms are there?
num_firms = size(temp_A,2);

lengths = zeros(num_firms*num_firms,1);

%Transpose A so that it is a normal adjacency matrix
A = temp_A';

%calculate the shortest path (technically min cost path) for A.
[cost,path] = dijkstra(A,A);


index = 1;
%for each firm
for i = 1:num_firms
    %to each firm
    for j = 1:num_firms
        %calculate the length of the shortest path
        path_length = length(path{i,j});
        %store the path length in the list of lengths
        lengths(index) = path_length;
        %increment index
        temp_index = index;
        index = temp_index + 1;
    end  
    
end

%get the maximum of the path lenghs
max_shortest_path = max(lengths);



end