%This function determines the average shortest path length for a network.
%It takes as input an adjacency matrix, Adj
%It returns the average shortest path length, avg_shortest_path

%This uses Joseph Kirk's "dijkstra.m" path finding algorithm.

function [avg_shortest_path] = avg_shortest_path_dist(Adj)

%Get the function variables
temp_A = Adj; %Ratchet adjacency matrix

%How many firms are there?
num_firms = size(temp_A,2);

avg_shortest_path = 0;

%Transpose A so that it is a normal adjacency matrix
A = temp_A';

sum_length = 0;

%calculate the shortest path (technically min cost path) for A.
[cost,path] = dijkstra(A,A);

%for each firm
for i = 1:num_firms
    %to each firm
    for j = 1:num_firms
        %calculate the length of the shortest path
        path_length = length(path{i,j});
        %add the length to the sum of the lengths
        temp_length = sum_length;
        sum_length = temp_length + path_length;
    end  
    
end

%divide by the number of paths = num_firms^2
avg_shortest_path = sum_length/(num_firms*num_firms);







end