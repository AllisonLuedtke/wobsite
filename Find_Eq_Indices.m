%this finds the indices of the eq networks

function [index_list] = Find_Eq_Indices(Eq_List, Feas_List)


%function variables
Eqs = Eq_List; %list of eq adjacency matrices
Feasibles = Feas_List; %list of feasible adjacency matrices

%how many firms
num_firms = size(Eqs,2);

%How many eq's?
num_eqs = length(Eqs)/num_firms;

index_list = zeros(num_eqs,1);

%for each eq, get the adj, find it in the list of feasibles
for e = 1:num_eqs
    stop_index = e*num_firms;
    start_index = stop_index - num_firms + 1;
    
    find_eq_mat = Eqs(start_index:stop_index,:);
    
    %find the index
    new_index = find_feasible_index(Feasibles,find_eq_mat);
    %store it
    index_list(e) = new_index;
    
    
    
end



end