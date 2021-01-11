%this simulation finds all of the pairwise stable and coordination proof
%equilibria, deletes each edge in each, and then finds them all again

num_firms = 5;
num_reps = 50;

%storage
num_cp_list = zeros(num_reps,1);
num_pw_list = zeros(num_reps,1);
ratio_list = zeros(num_reps,1);

for rep = 1:num_reps

    %initalize flag as 1
    flag = 1;
    %while the flag isn't switched to 0, keep drawing new A's and Z's until
    %you get one where the flag is switched to 0 (to ensure more than one
    %supplier for everyone)
    while (flag == 1)
        %adjacency matrix and productivity matrix
        [A, Z, f] = PCreate_A_and_Z(num_firms);
        flag = f;
    end
    
    %Construct the matrix of potential supplier indices
    sup_matrix = zeros(num_firms, num_firms - 1);
    
    for k = 1:num_firms
        %What are the indices of the nonzero elements (suppliers)?
        index_list = find(A(k,:));
        %How many are there?
        num_options = length(index_list);
        %For each of those options put the index in the supplier matrix
        for s = 1:num_options
            index = index_list(s);
            sup_matrix(k,s) = index;
        end
    end
    
    
    %find the original pwise and cpe solutions solutions
    [pwise_output, pwise_networks, pwise_res, x_star, Pay, eq_edge_list_pw, agg_output, org_eff, eff_list,...
        cp_output, cp_networks, cp_org_eff] = Model_Simulation_PW_and_CP_Mult(A,Z);
    %if the res is small enough and an equilibrium is found, set flag to 1
    if (pwise_res < 10^(-5)) && (pwise_networks(1,1) ~= 0) && (cp_networks(1,1) ~= 0)
       pwise_flag_A = 1;
    end
    
    %how many pw's are there?
    num_pw = size(pwise_networks,1);
    num_pw_list(rep) = num_pw;
    
    %how many cp's are there?
    num_cp = size(cp_networks,1);
    num_cp_list(rep) = num_cp;
    
    %how many of those pw's are cp?
    cp_ratio = num_cp/num_pw;
    ratio_list(rep) = cp_ratio;
   
    
    
end   