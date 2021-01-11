function [ F ] = SystemToSolve(V, Adj, prod, alpha, epsilon, Labor)
%This function creates the system of equations that determine a pairwise
%stable equilibrium

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

g = 0.5;

%create the matrix of choice variables from V
size_V = size(V);
last_index = size_V(1);
J = (last_index - 1)/6;  %number of firms

vector_of_x = V(1:last_index-1);
M = zeros(J,6);

%This should work but I don't have access to it.
%M = vec2mat(vector_of_x,6);
%THIS IS LIKE CRAZY INEFFICIENT
v_index = 1;
for h = 1:J
    for m = 1:6
        M(h,m) = vector_of_x(v_index);
        v_index = v_index + 1;
    end
end




w= V(last_index);

%The matrix X contains the choice variables for each firm:
%       [p_0(1)   y_0(1)   Mu(1)   Lam(1)   x_e(1)   l(1)]
%       [p_0(2)   y_0(2)   Mu(2)   Lam(1)   x_e(2)   l(2)]
%                       ...      ...
%       [p_0(J)   y_0(J)   Mu(J)   Lam(J)   x_e(J)   l(J)]

%Get the values from the inputs
%M = X;      %The matrix of firm specific choice variables
%w = wage;   %The wage


%The number of firms = the number of rows of X
%For now it's 3
%J = 3;

%last_index = J*6 + 1;
F = zeros(last_index,1);


%create the sum of p_j_0 ^ (1 - epsilon)
sum_p_j_e = 0;
for k = 1:J
    temp_p = M(k,1);
    sum_p_j_e = sum_p_j_e + (temp_p^(1 - e));
end

%Create the matrix of z's used in this particular network
%Eventually this will be the multiplication of the network being used
% and the big matrix of all the Z's
%Z = zeros(3);
%Z(1,2) = 10;
%Z(2,1) = 5;
%Z(3,1) = 1;

%Make a vector of productivities of edges pointing to each firm 
%(the z's for the specific network)
%Z(1) = 1/10;
%Z(2) = 1/5;
%Z(3) = 1;

%Adjacency Matrix

%For cycle of 3
%A(1,2) = 1;
%A(2,3) = 1;
%A(3,1) = 1;





%Calculate the prices for the particular network
%Eventually this will be a whole big thing.

%Vector of prices paid by each firm

P_E = Get_Prices(A,Z,w);


%P_E = Find_Prices_Two_Inputs(A, Z, a, g, w);

%Prices for cycle of 3
%P_E(3) = (1/(Z(1)))* (P_E(1)^a)*(w^(1 - a));
%P_E(2) = (1/(Z(3)))* (P_E(3)^a)*(w^(1 - a));

%P_E(1) = ((1/Z(2))*(1/(Z(1)^a))*(w^((a*(1-a)) + (1-a))))^(1/(1-a^2));
%P_E(2) = (1/(Z(1)))* (P_E(1)^a)*(w^(1 - a));
%P_E(3) = P_E(2);


%Create the equations for each firm
for j = 1:J
    i = j-1;
    
    %Set the variables so you're not using crazy indices
    p_e_j = P_E(j,1);
    z_e_j = Z(j);
    p_j_0 = M(j,1);
    y_j_0 = M(j,2);
    Mu_j = M(j,3);
    Lam_j = M(j,4);
    x_e_j = M(j,5);
    l_j = M(j,6);
    
    %get the vector of customers for firm j
    customers_j = get_customers(A,j);
    size_customers = size(customers_j);
    num_customers = size_customers(2);
    
    %get the sum of output used by customers of firm j
    sum_x_e = 0;
    if num_customers > 0
        for c = 1:num_customers
            cus = customers_j(c);
            sum_x_e = sum_x_e + M(cus,5);
        end
    end
    
    %set the final consumption demand for firm j
    %d_j = (w*L/ sum_p_j_e)*(p_j_0^(-e));
    
    %set the derivative of d_j
    %d_of_d_j = (w*L)*(((-e)*(p_j_0^(-e-1))*(sum_p_j_e) - ((1-e)*(p_j_0^(-2*e))))/(sum_p_j_e^2));
    
    %First Equation: p_j_0 FOC
    first_index = 1 + (i*6);
    %F(first_index) = y_j_0 + (Mu_j*d_of_d_j);
    F(first_index) = y_j_0*(sum_p_j_e*sum_p_j_e)*(p_j_0^(e+1)) + Mu_j*w*L*(((-e)*sum_p_j_e) - ((1 - e)*(p_j_0^(1 - e))));
    
    %Second Equation: y_j_0 FOC
    second_index = 2 + (i*6);
    F(second_index) = p_j_0 - Mu_j - Lam_j;
    
    %Third Equation: x_e_j FOC
    third_index = 3 + (i*6);
    %F(third_index) = -p_e_j + Lam_j*(a/((a^a)*((1 - a)^(1 - a))))*z_e_j*(x_e_j^(a-1))*(l_j^(1-a));
    F(third_index) = ((-p_e_j)*((a^a)*((1 - a)^(1 - a)))*(x_e_j^(1 - a))) + (a*Lam_j*z_e_j*(l_j^(1 - a)));
    
    %Fourth Equation: l_j FOC
    fourth_index = 4 + (i*6);
    %F(fourth_index) = 4;
    %F(fourth_index) = -w + Lam_j*((1 - a)/((a^a)*((1 - a)^(1 - a))))*z_e_j*(x_e_j^(a))*(l_j^(-a));
    F(fourth_index) = ((-w)*((a^a)*((1 - a)^(1 - a)))*(l_j^a)) + ((1 - a)*Lam_j*z_e_j*(x_e_j^a));
    
    
    %Fifth Equation: final consumption output constraint
    fifth_index = 5 + (i*6);
    %F(fifth_index) = y_j_0 - d_j;
    F(fifth_index) = (y_j_0*(p_j_0^e)*(sum_p_j_e)) - w*L;
    
    %Sixth Equation: total output constraint
    sixth_index = 6 + (i*6);
    %F(sixth_index) = y_j_0 + sum_x_e - (1/((a^a)*((1 - a)^(1 - a))))*z_e_j*(x_e_j^(a))*(l_j^(1-a));
    F(sixth_index) = (y_j_0 + sum_x_e)*((a^a)*((1 - a)^(1 - a))) - z_e_j*(x_e_j^(a))*(l_j^(1-a));
end


%The labor clearing constraint
%Sum across last column of M (=X)
matrix_sum = sum(M);
l_used = matrix_sum(6);

F(last_index) = l_used - L;     %Labor Clearing







end

