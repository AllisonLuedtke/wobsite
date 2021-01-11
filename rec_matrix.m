%This is a recurseive function that enumerates all of the different
%   combinations of suppliers possible for a given potential network
%   adjacency matrix. 

%It takes as input
%   - firm (row) we're on
%   - a working vector naming the suppliers used
%   - the potential adjacency matrix
%   - matrix to store the finished list for each possible eq. network

function [matrix] = rec_matrix(f, sup_vec, sup_mat, mat)

firm = f; %which firm are we on
supplier_vector = sup_vec; %list of suppliers we're building
supplier_matrix = sup_mat; %Matrix describing the suppliers available to 
                            %each firm
matrix = mat;   %matrix where we're storing the final product

%How many firms are we working with?
num_firms = size(matrix,2);



%RETURN CONDITION

%If we're on the last firm
if firm == num_firms
    
    %Print the list of supplier's we've created on this round
    %supplier_vector
    
    %What's the next open spot in the matrix?
    next_spot = find(matrix(:,1)==0,1,'first'); 
    
    %Store it in the matrix
    %Only make changes to this matrix if you're on the last firm and have
    %completed a list of suppliers
    matrix(next_spot,:) = supplier_vector;
    
    %Return
    
    
else
%Otherwise, for each potential supplier of the next firm, call the function


    %What are the supplier options for the next firm?
    firm_suppliers = nonzeros(supplier_matrix(firm+1,:));
    %How many?
    n = length(firm_suppliers);

    %For each such supplier, add him to the vector of suppliers and call the
    %function with it.
    for i = 1:n
        %create a temp vector so you don't ruin it
        temp_vec = supplier_vector;
        %put potential supplier i in the supplier vector
        temp_vec(firm+1) = firm_suppliers(i);
        
        %call the function with the new working vector for the next firm
        matrix = rec_matrix(firm+1,temp_vec,supplier_matrix,matrix);
    
    end


end

