%this function finds the index of a given adjacency matrix for a feasible
%network among a list of all the possible feasible networks

function [index] = find_feasible_index(f_list, find_feas)

%get the function variables
feasible_list = f_list;  %the list of all possible feasible networks
other = find_feas;  %the network we want to find in the list

%how many firms?
num_firms = size(other,1);

%how many feasible networks?
num_feas = size(feasible_list,1)/num_firms;

%initialize index
index = -1;

for feas = 1:num_feas
    n = feas - 1;
    start = n*num_firms + 1;
    stop = start + num_firms - 1;
    
    adj = feasible_list(start:stop,:);
    
    %if the adj is equal to the feasible network we're checking this is its
    %index
    
    if prod(prod(adj == other)) == 1
        %they're equal
        index = feas;
    end
  
end

end