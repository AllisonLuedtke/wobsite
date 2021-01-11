function [ F ] = Planner_System(V, Adj, prod, alpha, epsilon, Labor)
%This function creates the system of equations that determine the solution
%to the Planner's Problem

%It takes a possible equilibrium network adjacency matrix

%It takes as inputs:
%   - V, the list of choice variables, turned into a matrix here
%   - Adj, the adjacency matrix of the network in question,
%   - prod, the list of productivities for the edge pointing to each firm 
%   - alpha
%   - epsilon
%   - L

%Get the function values
A = Adj;
Z = prod;
a = alpha;
e = epsilon;
L = Labor;

%create the matrix of choice variables from V
size_V = size(V);
last_index = size_V(1);
J = (last_index - 1)/4;  %number of firms

vector_of_x = V(1:last_index-1);
M = zeros(J,4);

%This should work but I don't have access to it.
%M = vec2mat(vector_of_x,6);
%THIS IS LIKE CRAZY INEFFICIENT
v_index = 1;
for h = 1:J
    for m = 1:4
        M(h,m) = vector_of_x(v_index);
        v_index = v_index + 1;
    end
end
%M



mu= V(last_index);

%The matrix M contains the choice variables for each firm:
%       [y_0(1)   Lam(1)   x_e(1)   l(1)]
%       [y_0(2)   Lam(1)   x_e(2)   l(2)]
%                       ...      ...
%       [y_0(J)   Lam(J)   x_e(J)   l(J)]


F = zeros(last_index,1);

%Create the function for each firm
%This will fill the values of F from 1 to J-1

for j = 1:J
    i = j-1;
    
    %Set the variables so you're not using crazy indices
    y_j_0 = M(j,1);
    lam_j = M(j,2);
    x_e_j = M(j,3);
    l_j = M(j,4);
    z_e_j = Z(j);
    
    %What is the index of the supplier of j?
    supplier_j = find(A(j,:));
    lam_s = M(supplier_j,2);
    
    %get the vector of customers for firm j
    customers_j = get_customers(A,j);
    size_customers = size(customers_j);
    num_customers = size_customers(2);
    
    %get the sum of output used by customers of firm j
    sum_x_e = 0;
    if num_customers > 0
        for c = 1:num_customers
            cus = customers_j(c);
            temp_x_e = sum_x_e;
            sum_x_e = temp_x_e + M(cus,3);
        end
    end
    
    
    %First Equation: y_j_0 FOC
    first_index = 1 + (i*4);
    %create sum of y_j_o ^(e-1/e)
    sum_y_k_e = sum(M(:,1).^((e-1)/e));
    F(first_index) = 1 - (lam_j*(y_j_0 ^ (1/e))*(sum_y_k_e ^ (1/(1-e))));
    
    %Second Equation: l_j FOC
    second_index = 2 + (i*4);
    F(second_index) = (1-a)*lam_j*z_e_j*(x_e_j^a) - mu*((a^a)*((1-a)^(1-a)))*(l_j^a);
    
    %Third Equation: x_e_j FOC
    third_index = 3 + (i*4);
    F(third_index) = a*lam_j*z_e_j*(l_j^(1-a)) - lam_s*((a^a)*((1-a)^(1-a)))*(x_e_j^(1-a));
    
    %Fourth Equation: production constraint
    fourth_index = 4 + (i*4);
    F(fourth_index) = (((a^a)*((1-a)^(1-a)))*(y_j_0 + sum_x_e)) - (z_e_j*(x_e_j^a)*(l_j^(1-a)));
    
end

%Labor clearing constraint
sum_l_j = sum(M(:,4));
F(last_index) = L - sum_l_j;

end




