function [Eq_Network, residuals, num_mats] = Model_Simulation_Mult_Eq(Adj,Prod)
%This runs the simulation of my model (finds the pairwise stable eq)

%Get the Potential Network Adjacancy matrix and the productivity matrix.
%Eventually you will need to write code to CREATE these two things. But for
%now I'm just hardcoding it.\

%Get function variables
A = Adj;
Z = Prod;

%A = [0 1 1; 1 0 1; 0 1 0]
%Original
%A = [0 1 0 1; 1 0 1 0; 1 1 0 0; 1 0 1 0]
%Variations
%A = [0 0 0 1; 1 0 1 0; 1 1 0 0; 1 0 1 0]
%A = [0 1 0 1; 0 0 1 0; 1 1 0 0; 1 0 1 0]
%A = [0 1 0 1; 1 0 1 0; 0 1 0 0; 1 0 1 0]
%A = [0 1 0 1; 1 0 1 0; 1 1 0 0; 1 0 0 0]
%Z = [0 3 0 2; 5 0 6 0; 11 10 0 0; 8 0 7 0]
%Z = [0 5 0 2; 10 0 12 0; 20 19 0 0; 29 0 30 0]
%Z = ones(4,4)
%Z = [0.2316 0.3955 0.8852 0.2619; 0.4889 0.3674 0.9133 0.3354; 0.6241 0.9880 0.7962 0.6797; 0.6791 0.0377 0.0987 0.1366]
%Z = [2 3 1 6; 3 9 5 5; 8 8 5 2; 2 6 1 9]
%Z = [0 1/10 1/10; 1/5 0 1/5; 1 1 0]

%Set the model scalars
a = 0.33;
e = 0.10;
L = 1;

%How many firms
num_firms = size(A,1);

%Get all of the possible equilbrium networks for that potential network
Eq_Mat_List = Find_Potentials(A);
num_mats = size(Eq_Mat_List,1);

%Create an empty matrix for the payoffs
Payoffs = zeros(num_mats,num_firms);

Aggregate_Outputs = zeros(num_mats,1);
Other_Agg_Outputs = zeros(num_mats,1);
Efficiencies = zeros(num_mats,num_firms);


%For each possible equilibrium, make it into an adjacency matrix and z's
for m = 1:num_mats
    %turn the row of the matrix into an adjacency matrix
    %create a blank matrix
    Eq_A = zeros(num_firms);
    
    %get the row
    row = Eq_Mat_List(m,:);
    
    %for each firm, put a 1 in the place of it's supplier in the adjacency
    %matrix
    for f = 1:num_firms
        customer = f;
        supplier = row(f);
        Eq_A(customer,supplier) = 1;
    end
    
    %construct the vector of z's
    %by multiplying Eq_A by the Z matrix (element by element) then a column of ones
    Eq_z = (Eq_A.*Z)*ones(num_firms,1);
    
    %Get the edge prices
    %edge_prices = Get_Prices(Eq_A,Eq_z);
    
    %Find Firm Choices
    %fun = @SystemToSolve;
    guess_size = num_firms*6 + 1;
    guess = 0.1*ones(guess_size,1);
    lower = zeros(guess_size,1);
    %for b = 1:num_firms
     %   a = b - 1;
      %  bound_index_1 = a*6 + 3;
       % bound_index_2 = a*6 + 4;
        %lower(bound_index_1) = -1000000000000;
        %lower(bound_index_2) = -1000000000000;  
    %end
  
    
    [lsq_guess, res] = lsqnonlin(@(v) SystemToSolve(v,Eq_A,Eq_z,a,e,L),guess,lower);
    %[x_star, fval] = fsolve(@(v) SystemToSolve(v,Eq_A,Eq_z,a,e,L),lsq_guess)
    x_star = lsq_guess;
    residuals = res;
    
    %construct the matrix of firm choices from the vector x_star
    X_Mat = zeros(num_firms,6);
    x_index = 1;
    for h = 1:num_firms
        for g = 1:6
            X_Mat(h,g) = x_star(x_index);
            temp_index = x_index;
            x_index = temp_index + 1;
        end
    end
    %X_Mat
    
    w_star = x_star(end);

    
    
    %construct payoffs for each firm
    
    %get the prices
    
    %******THERE"S PROBABLY A BETTER WAY TO DO THIS*******
    
    edge_price = Get_Prices(Eq_A,Eq_z,w_star);
    %Payoffs

    for j = 1:num_firms
        
        payoff = 0;
        
        %cons_rev = p_j0 * y_j0
        cons_rev = X_Mat(j,1)*X_Mat(j,2);
        %intermediate input revenue
        %find the customers for firm j
        cus = get_customers(Eq_A,j);
        num_cus = length(cus);
        network_rev = 0;
        
        %for each customer, add the revenue from them
        for c = 1:num_cus
            %network_rev = network_rev + p(e_c)*x(e_c)
            temp_rev = network_rev;
            network_rev = temp_rev + edge_price(c)*X_Mat(c,5);
        end
        
        %calculate the cost of firm j's intermediate input
        %network_cost = p(e_j)*x(e_j)
        network_cost = edge_price(j)*X_Mat(j,5);
        
        %labor_cost = w_star*l_j
        labor_cost = w_star*X_Mat(j,6);
        
        payoff = cons_rev + network_rev - network_cost - labor_cost;

        Payoffs(m,j) = payoff;
        
        
    end
    
    %calculate aggregate output for that equilibrium matrix
    
    %first create a vector of efficiencies
    q_j = Eq_z.*(edge_price.^(-a))*(w_star^a);
    Efficiencies(m,:) = q_j;
    
    %create the weighted q's
    q_j_e = q_j.^(e - 1);
    sum_q_j_e = sum(q_j_e);
    
    %calculate aggregate output
    Agg_Output = (sum_q_j_e^(1/(e - 1)))*L;
    Aggregate_Outputs(m) = Agg_Output;
    
    %calculate the other measure of aggregate output
    y_j_e = X_Mat(:,2).^((e-1)/e);
    Other_Agg_Output = (sum(y_j_e))^(e/(e-1));
    Other_Agg_Outputs(m) = Other_Agg_Output;
    
    
    
end
%Payoffs

%Construct the matrix of potential supplier indices
sup_matrix = zeros(num_firms, num_firms - 1);
for i = 1:num_firms
    %What are the indices of the nonzero elements (suppliers)?
    index_list = find(A(i,:));
    %How many are there?
    num_options = length(index_list);
    %For each of those options put the index in the supplier matrix
    for s = 1:num_options
        index = index_list(s);
        sup_matrix(i,s) = index;
    end    
end
%sup_matrix

%Pick the equilibrium
%Eq_Mat_List
%Payoffs

%try rounding the payoffs
%Temp_Payoffs = round(10000*Payoffs)/10000;

%initalize Eq_Output and Also Output
Eq_Output = -1;
Also_Output = -1;

%sup_matrix
Eq_Network = Pick_Equilibrium_Mult(Eq_Mat_List,Payoffs,sup_matrix);




end

