%This is a recursive function that visits each node of a branch and
%calculates is price. It takes as inputs:
%       - a root node index, r
%       - a list of customers, cus
%       - An adjacency matrix, Adj
%       - a productivity matrix, Prod
%       - A price vector, p_e
%       - alpha, the alpha
%       - wage, the wage

%It returns the new and updated price vector, p_e_new

function [p_e_new] = rec_branch_price(r, cus, Adj, Prod, p_e, alpha, wage)

%Get the function values
root = r; %index of the root node
customers = cus; %list of customers
A = Adj; %Adjacency matrix
z = Prod; % vector of productivities

p_e_temp_1 = p_e; %vector of prices to be built throughout once you add the new stuff

a = alpha; %alpha
w = wage; %wage

if isempty(customers) %If no customers, return
    p_e_new = p_e_temp_1;
    return;   
else
    %calculate the MC of the root
    MC = (1/z(root))*(p_e_temp_1(root)^a)*(w^(1 - a));
    
    num_customers = length(customers);
    
    %for each customer:
    for c = 1:num_customers
        cus_index = customers(c);
        p_e_temp_1(cus_index) = MC;
        
        %find the customers of the current customer
        next_customers = find(A(:,cus_index));
        
        %call function for the customer
        p_e_temp_2 = rec_branch_price(cus_index,next_customers,A,z,p_e_temp_1,a,w);
    
    end
    
    
end
p_e_new = p_e_temp_2;
end