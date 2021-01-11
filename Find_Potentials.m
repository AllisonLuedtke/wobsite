%This function takes as input a potential network adjacency matrix, P, and
%returns a matrix, Pot_Eq, containing all the possible lists of suppliers (that
%would make up an equilibrium network's adjacency matrix). It calls the
%recursive function rec_matrix.

function [Pot_Eq] = Find_Potentials(P)

%get function variables
Potential = P;

%How many firms?
num_firms = size(Potential,1);

%Construct the matrix of potential supplier indices
sup_matrix = zeros(num_firms, num_firms - 1);
%And count how many rows you'll need for the storing matrix.
num_rows = 1;
for i = 1:num_firms
    %What are the indices of the nonzero elements (suppliers)?
    index_list = find(P(i,:));
    %How many are there?
    num_options = length(index_list);
    temp_rows = num_rows;
    num_rows = temp_rows*num_options;
    %For each of those options put the index in the supplier matrix
    for s = 1:num_options
        index = index_list(s);
        sup_matrix(i,s) = index;
    end    
end
%sup_matrix

%construct the temporary storing matrix
temp_store = zeros(num_rows,num_firms);

%construct initial vector
init_vector = zeros(num_firms,1);

%call the function for firm 0
return_mat = rec_matrix(0,init_vector,sup_matrix,temp_store);

Pot_Eq = return_mat;

end




