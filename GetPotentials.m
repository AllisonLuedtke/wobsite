%This function takes as input a potential network adjacency matrix (P_A)
% and returns a list of of all possible equilibrium networks (Pot_A)

function [  ] = GetPotentials(P_A)

Potential = P_A; %Potential Adjacency Matrix
num_firms = size(Potential,1);

%create a vector of the num of suppliers for each firm
num_suppliers = zeros(num_firms,1);

for i = 1:num_firms
    num_suppliers(i) = nnz(Potential(i,:));
end
num_suppliers

%Find the total number of potential matrices(= product of num suppliers)
%                                       = prod of num_suppliers entries

num_potentials = prod(num_suppliers);
num_rows = num_potentials*num_firms;

Big_Matrix = zeros(num_rows,num_firms);

%For each firm, set the row of each of the potential matrices

for r = 1:num_firms
    %Set half or a third or whatever to be each of the potential suppliers
    s = num_suppliers(r);
    
    
end



