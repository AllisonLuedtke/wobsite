%This script runs the simulation that compares the planner's solution to
%pairwise stable and the before and after output when you delete an edge

%Firms never have fewer than two suppliers in this version

%set the number of repetitions
num_reps = 5;
%set number of firms
num_firms = 5;

%initialize storage vectors (for each rep) (each size num_reps x 1)
planner_flag_A = zeros(num_reps,1); %flag indicating whether you found a planner solution
pwise_flag_A = zeros(num_reps,1);   %flag indicating whether you find a pwise solution ***NOTE that the pwise flag applies for cp too
connectivity_A = zeros(num_reps,1); %avg shortest path distance
planner_connectivity_E = zeros(num_reps,1); %avg shortest path distance of the planner solution eq
pwise_connectivity_E = zeros(num_reps,1); %avg hsortest path distance of the pwise eq
cp_connectivity_E = zeros(num_reps,1); %avg shortest path distance of the original cp eq
planner_relative_output = zeros(num_reps,1); %pwise output/planner output
cp_planner_relative_output = zeros(num_reps,1); %cp output/planner output
cp_pwise_relative_output = zeros(num_reps,1); %pwise output/ cp output
avg_planner_relative = zeros(num_reps,1); %average relative output for the planner and edge deletion
avg_pwise_relative = zeros(num_reps,1); %average relative output for pwise and edge deletion

%%initialize storage vectors (for each edge deletion) (each size
%%num_reps*num_firms x 1)
num_edges_deleted = num_reps*num_firms;
planner_relative = zeros(num_edges_deleted,1); %relative output for the planner and edge deletion
new_output_pw = zeros(num_edges_deleted,1); %holds the new output after the edge deletion pw
new_output_pl = zeros(num_edges_deleted,1); %holds the new output after the edge deletion pl
new_output_cp = zeros(num_edges_deleted,1); %holds the new output after the edge deletion cp
new_num_cust_pl = zeros(num_edges_deleted,1); %holds the number of customers of the customer of the deleted edge after
new_num_cust_pw = zeros(num_edges_deleted,1); %for pwise
new_num_cust_cp = zeros(num_edges_deleted,1); %for cp
new_pl_conn_E = zeros(num_edges_deleted,1); %holds the connectivity for the new network for the pl
new_pw_conn_E = zeros(num_edges_deleted,1); %holds the connectivity for the new network for the pw
new_cp_conn_E = zeros(num_edges_deleted,1); %holds the connectivity for the new network for the cp
pwise_relative = zeros(num_edges_deleted,1);   %relative output for pwise and edge deletion
cp_relative = zeros(num_edges_deleted,1);   %relative output for cp edge deletion
planner_poss_sup = zeros(num_edges_deleted,1); %number of possible suppliers for each customer of deleted edges in the planner experiment
pwise_poss_sup = zeros(num_edges_deleted,1); %num poss suppliers for each customer in the pwise experiment
cp_poss_sup = zeros(num_edges_deleted,1);   %num possible suppliers for each customer in the cp experiment
planner_num_cust = zeros(num_edges_deleted,1); % number of customers in eq of the customer of the deleted edge in the planner experiemnt
pwise_num_cust = zeros(num_edges_deleted,1); %num customers of the customer in the pwise experiment
cp_num_cust = zeros(num_edges_deleted,1);   %num customers of the customer in the cp experiment
planner_flag_edge = zeros(num_edges_deleted,1); %flag indicating whether you found a planner solution in the experiment
pwise_flag_edge = zeros(num_edges_deleted,1);   %flag indicating whether you find a pwise solution in the experiment
cp_flag_edge = zeros(num_edges_deleted,1);      %flag indicating whether you find a cp solution in the experiment
big_planner_connectivity_E = zeros(num_edges_deleted,1); %vector containing the connectivity of the planner eq network at the edge deletion level
big_pwise_connectivity_E = zeros(num_edges_deleted,1); %vector containing the connectivity of the pwise eq network at the edge deletion level
big_cp_connectivity_E = zeros(num_edges_deleted,1); %vector containing the connectivity of the cp eq network at the edge deletion level
big_planner_connectivity_E_new = zeros(num_edges_deleted,1); %vector containing the connectivity of the planner eq network at the edge deletion level
big_pwise_connectivity_E_new = zeros(num_edges_deleted,1); %vector containing the connectivity of the pwise eq network at the edge deletion level
big_cp_connectivity_E_new = zeros(num_edges_deleted,1); %vector containing the connectivity of the cp eq network at the edge deletion level

%vectors for who bares the cost
original_pwise_eff = zeros(num_edges_deleted, num_firms); %holds the list of efficiencies for the original pwise solution
new_pwise_eff = zeros(num_edges_deleted, num_firms); %holds the list of efficiencies for the new pwise solution
j_star_old_eff_pw = zeros(num_edges_deleted, 1); %hold the efficiency of j* before
j_star_new_eff_pw = zeros(num_edges_deleted, 1); %hold the efficiency of j* after
cust_old_eff_pw = zeros(num_edges_deleted, num_firms); %holds the efficiencies of the customers of j* before 
cust_new_eff_pw = zeros(num_edges_deleted, num_firms); %holds the efficiencies of the customers of j* after

%counters for storing the edge deletion values
pl_index = 1;
pw_index = 1;


%big loop
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
    
    %find average shortest path distance
    avg_dist = avg_shortest_path_dist(A);
    %store it
    connectivity_A(rep) = avg_dist;
    
    %Calculate number of potential suppleirs for each firm
    num_pot_sups = sum(A~=0,2);
    

   
    
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
    
    %find the original pwise solution
    [pwise_outputs, pwise_networks, pwise_res, x_star, Pay, eq_edge_list_pw, agg_output, org_eff, eff_list,...
        cp_outputs, cp_networks, cp_org_eff] = Model_Simulation_PW_and_CP_Mult(A,Z);
    %if the two networks are connected components
    num_pw_cycs = 0;
    if pwise_networks(1,1) ~= 0
        %pick the one with the highest output to be THE one
        [max_pw_out, max_pw_out_i] = max(pwise_outputs);
        pwise_network = pwise_networks(max_pw_out_i,:);
        pwise_output = max_pw_out;
        
        
    %turn each into an Adj matrix, check cycles
    pwise_adj = zeros(num_firms);
    %for each firm, put a 1 in the place of it's supplier in the adjacency
    %matrix
    for fr = 1:num_firms
        customer_pw = fr;
        supplier_pw = pwise_network(fr);
        pwise_adj(customer_pw,supplier_pw) = 1;
    end
    num_pw_cycs = size(Find_Cycle_For_Test(pwise_adj),1);
    end
    %first turn the cp eq network
    num_cp_cycs = 0;
    if cp_networks(1,1) ~= 0
        %pick the one with the highest output
        [max_cp_out, max_cp_out_i] = max(cp_outputs);
        cp_network = cp_networks(max_cp_out_i,:);
        cp_output = max_cp_out;
        
    cp_adj = zeros(num_firms);
    %for each firm put a 1 in for its supplier
    for cp_fr = 1:num_firms
        customer_cp = cp_fr;
        supplier_cp = cp_network(cp_fr);
        cp_adj(customer_cp,supplier_cp) = 1;
    end
    num_cp_cycs = size(Find_Cycle_For_Test(cp_adj),1);
    end
    
    %if the res is small enough and an equilibrium is found and its connected, set flag to 1
    if (pwise_res < 10^(-5)) && (pwise_network(1,1) ~= 0) && (cp_network(1,1) ~= 0) && (num_pw_cycs == 1) && (num_cp_cycs == 1)
       pwise_flag_A(rep) = 1; 
    end

    
    %IF FLAG DO PLANNER
    if pwise_flag_A(rep) == 1
    
    %find the avg shortest path distance for the pwise network
    
    %find the average shortest path for that adjacency matrix
    avg_dist_pw = avg_shortest_path_dist(pwise_adj);
    %store it
    pwise_connectivity_E(rep) = avg_dist_pw;
    
    %Find the number of customers used in the original pwise equilibrium
    num_pwise_cust = zeros(num_firms,1);
    for c = 1:num_firms
        %How many times does firm c appear as a supplier in the eq network
        num_cus = sum(pwise_network == c);
        num_pwise_cust(c) = num_cus;
    end
    
    %%%FOR CP%%%%%
    
    %find the average shortest path distance for that adj
    avg_dist_cp = avg_shortest_path_dist(cp_adj);
    %store it
    cp_connectivity_E(rep) = avg_dist_cp;
    
    %find the number of customers for each firm in this original cp eq
    num_cp_cust = zeros(num_firms,1);
    for cp_c = 1:num_firms
        num_cus_cp = sum(cp_network == cp_c);
        num_cp_cust(cp_c) = num_cus_cp;
    end
    
    
    
    %construct the planner guess from the p_wise solution
    planner_guess = zeros(num_firms*4 + 1,1);
    for f = 1:num_firms
       i = f-1;
       pwise_index = i*6;
       planner_index = i*4;
       %y_j_0 guess
       planner_guess(planner_index + 1) = x_star(pwise_index + 2);
       %lam_j guess
       planner_guess(planner_index + 2) = x_star(pwise_index + 4);
       %x_e_j guess
       planner_guess(planner_index + 3) = x_star(pwise_index + 5);
       %l_j guess
       planner_guess(planner_index + 4) = x_star(pwise_index + 6);    
    end
    
    %find the original planner solution
    [planner_output, planner_network, planner_res, Ut, eq_edge_list_pl] = Planner_Simulation(A,Z, planner_guess);
    %check connected
    %first turn it into an adjacency matrix
    planner_adj = zeros(num_firms);
    %for each firm, put a 1 in the place of it's supplier in the adjacency
    %matrix
    for frm = 1:num_firms
        customer_pl = frm;
        supplier_pl = planner_network(frm);
        planner_adj(customer_pl,supplier_pl) = 1;
    end
    num_pl_cycs = size(Find_Cycle_For_Test(planner_adj),1);
    
    if (planner_res < 10^(-5)) && (num_pl_cycs == 1)
       planner_flag_A(rep) = 1; 
    end
    
    %find the avg shortest path distance for the planner network
    %first turn it into an adjacency matrix
    planner_adj = zeros(num_firms);
    %for each firm, put a 1 in the place of it's supplier in the adjacency
    %matrix
    for frm = 1:num_firms
        customer_pl = frm;
        supplier_pl = planner_network(frm);
        planner_adj(customer_pl,supplier_pl) = 1;
    end
    
    %find the avg shortest path dist
    avg_dist_pl = avg_shortest_path_dist(planner_adj);
    %store it
    planner_connectivity_E(rep) = avg_dist_pl;
    
    
    %Find the number of customers used in the original planner equilibrium
    num_planner_cust = zeros(num_firms,1);
    for c = 1:num_firms
        %How manh times does firm c appear as a supplier in the eq network
        num_cus = sum(planner_network == c);
        num_planner_cust(c) = num_cus;
    end
    
    %END FOR IF PWISE SUCCESSFUL
    end
    
    %if both planner and pwise successful, then bother with the rest
    
    %IF BOTH FLAGS == 1
    if (pwise_flag_A(rep) == 1) && (planner_flag_A(rep) == 1)
    
    %for each edge in the pwise solution, delete it calculate new eq
    for pw_c = 1:num_firms
        %get the supplier of the pw_c'th firm in the original equilibrium
        pw_edge_supplier = pwise_network(pw_c);
        
        %record that edge's customer's possible suppliers and number of
        %customers in eq
        pwise_poss_sup(pw_index) = num_pot_sups(pw_c);
        pwise_num_cust(pw_index) = num_pwise_cust(pw_c);
        
        %delete the payoffs, agg output, and eq_edge_lists that correspond to that
        %particular edge
        
        %anywhere a pw_edge_supplier appears in the pw_c column, delete
        %that row in them all.
        %this is the list of rows to delete
        delete_index_list = find(eq_edge_list_pw(:,pw_c)==pw_edge_supplier);
        
        %for each index, remove that row from the eq_edge_list, the agg output list
        %and from the payoff matrix
        
        %first copy them
        new_edge_list_pw = eq_edge_list_pw;
        new_payoffs = Pay;
        new_agg_output = agg_output;
        new_eff_list = eff_list;
        
        %delete the rows that correspond to the delete list
        new_edge_list_pw(delete_index_list,:) = [];
        new_payoffs(delete_index_list,:) = [];
        new_agg_output(delete_index_list,:) = [];
        new_eff_list(delete_index_list,:) = [];
        
        
        %create new sup list
        %copy A
        new_A = A;
        new_A(pw_c,pw_edge_supplier) = 0;
        %redo the sup list creating process with this new A
        new_sup_matrix = zeros(num_firms, num_firms - 1);
        for u = 1:num_firms
            %What are the indices of the nonzero elements (suppliers)?
            new_index_list = find(new_A(u,:));
            %How many are there?
            new_num_options = length(new_index_list);
            %For each of those options put the index in the supplier matrix
            for t = 1:new_num_options
                new_index = new_index_list(t);
                new_sup_matrix(u,t) = new_index;
            end    
        end
        
        %find the new eq
        new_eq_pws = Pick_Equilibrium_Mult(new_edge_list_pw,new_payoffs,new_sup_matrix);
        
        %if you successfully find the eq, do this and set the flag to 1
        if new_eq_pws(1,1) ~= 0
            %set the one with the largest output to be THE one
            %find the indices
            new_eq_indices = find(ismember(new_edge_list_pw,new_eq_pws,'rows'));
            %get the associated outputs
            new_eq_outputs = new_agg_output(new_eq_indices);
            %take the max
            [max_new_eq_pw, max_new_eq_pw_i] = max(new_eq_outputs);
            new_eq_pw = new_eq_pws(max_new_eq_pw_i,:);
            
        %find the new number of customers in the new eq network
        %of the deleted edge
        %How many times does the new edge appear in the new eq network?
        new_num_cus = sum(new_eq_pw == pw_c);
        %store it
        new_num_cust_pw(pw_index) = new_num_cus;
        
        %find the connectivity of the new network
        %first, turn it into an adjacency matrix
        new_pwise_adj = zeros(num_firms);
        %for each firm, put a 1 in the place of it's supplier in the adjacency
        %matrix
        for frm_pw = 1:num_firms
            new_cust_pw = frm_pw;
            new_sup_pw = new_eq_pw(frm_pw);
            new_pwise_adj(new_cust_pw,new_sup_pw) = 1;
        end
        %find avg shortest path distance
        new_pw_path_dist = avg_shortest_path_dist(new_pwise_adj);
        %store it
        new_pw_conn_E(pw_index) = new_pw_path_dist;
        
            
        %find the index of that new eq network
        new_eq_index = find(ismember(new_edge_list_pw,new_eq_pw,'rows'));
        
        %find the new agg output measure
        new_pwise_output = new_agg_output(new_eq_index);
        
        %store the new output
        new_output_pw(pw_index) = new_pwise_output;
        
        %calculate this relative output and store it
        pw_relative_output = new_pwise_output/pwise_output;
        pwise_relative(pw_index) = pw_relative_output;
        
        %store the connectivity
        big_pwise_connectivity_E(pw_index) = avg_dist_pw;
        
        %find the associate new list of efficiencies
        new_effs = new_eff_list(new_eq_index);
        %store the original effficiencies
        original_pwise_eff(pw_index,:) = org_eff;
        %store the new efficiencies
        new_pwise_eff(pw_index,:) = new_effs;
        
        %store the old efficiency of j*
        old_j_star_eff = original_pwise_eff(pw_index,pw_c);
        j_star_old_eff_pw(pw_index) = old_j_star_eff;
        
        %store the new eff of j*
        new_j_star_eff = new_pwise_eff(pw_index, pw_c);
        j_star_new_eff_pw(pw_index) = new_j_star_eff;
        
        %store the old and new effs of the customers of j*
        cust_index = 1;
        for p = 1:num_firms
            %if pw_c is the supplier of firm p in the original network
            if (pw_c == pwise_network(p))
                %store p's old efficiency
                %pw_index
                %cust_index
                old_cust_q = original_pwise_eff(pw_index,p);
                cust_old_eff_pw(pw_index,cust_index) = old_cust_q;
                %store p's new efficiency
                new_cust_q = new_pwise_eff(pw_index, p);
                cust_new_eff_pw(pw_index,cust_index) = new_cust_q;
                %increment the customer index
                temp_cust_index = cust_index;
                cust_index = temp_cust_index + 1;
            end   
        end
        
        
        
        %set big flag to 1
        pwise_flag_edge(pw_index) = 1;
        
        %SUCCESS FLAG END
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now do the Coalition-Proof Stuff!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Note that I'm using the pw_c counter
        cp_edge_supplier = cp_network(pw_c);
        
        %record that edge's customer's possible suppliers and number of
        %customers in eq
        cp_poss_sup(pw_index) = num_pot_sups(pw_c);
        cp_num_cust(pw_index) = num_cp_cust(pw_c);
        
        %anywhere a cp_edge_supplier appears in the pw_c column, delete
        %that row in them all.
        %this is the list of rows to delete
        cp_delete_index_list = find(eq_edge_list_pw(:,pw_c)==cp_edge_supplier);
        
        
        %for each index, remove that row from the eq_edge_list, the agg output list
        %and from the payoff matrix
        
        %first copy them
        new_edge_list_cp = eq_edge_list_pw;
        new_payoffs_cp = Pay;
        new_agg_output_cp = agg_output;
        new_eff_list_cp = eff_list;
        
        %delete the rows that correspond to the delete list
        new_edge_list_cp(cp_delete_index_list,:) = [];
        new_payoffs_cp(cp_delete_index_list,:) = [];
        new_agg_output_cp(cp_delete_index_list,:) = [];
        new_eff_list_cp(cp_delete_index_list,:) = [];
        
        %create new sup list
        %copy A
        new_cp_A = A;
        new_cp_A(pw_c,cp_edge_supplier) = 0;
        %redo the sup list creating process with this new A
        new_cp_sup_matrix = zeros(num_firms, num_firms - 1);
        for cp_u = 1:num_firms
            %What are the indices of the nonzero elements (suppliers)?
            cp_new_index_list = find(new_cp_A(cp_u,:));
            %How many are there?
            cp_new_num_options = length(cp_new_index_list);
            %For each of those options put the index in the supplier matrix
            for cp_t = 1:cp_new_num_options
                cp_new_index = cp_new_index_list(cp_t);
                new_cp_sup_matrix(cp_u,cp_t) = cp_new_index;
            end    
        end
        
        %find the new cp eq
        new_eq_cps = Find_CP_Eq_Mult(new_edge_list_cp,new_payoffs_cp,new_cp_sup_matrix);
        
        %if you successfully find the eq, do this and set the flag to 1
        if new_eq_cps(1,1) ~= 0
            %set the one with the max output to be the one
            %find the indices
            new_cp_eq_indices = find(ismember(new_edge_list_cp,new_eq_cps,'rows'));
            %get the outputs
            new_eq_cp_outputs = new_agg_output_cp(new_cp_eq_indices);
            %find the max
            [max_new_cp_eq_out, max_new_cp_eq_out_i] = max(new_eq_cp_outputs);
            %set it
            new_eq_cp = new_eq_cps(max_new_cp_eq_out_i,:);
            
            
            
        %find the new number of customers in the new cp eq network
        %of the deleted edge
        %How many times does the new edge appear in the new eq network?
        new_cp_num_cus = sum(new_eq_cp == pw_c);
        %store it
        new_num_cust_cp(pw_index) = new_cp_num_cus;
        
        %find the connectivity of the new network
        %first, turn it into an adjacency matrix
        new_cp_adj = zeros(num_firms);
        %for each firm, put a 1 in the place of it's supplier in the adjacency
        %matrix
        for cp_firm = 1:num_firms
            new_cust_cp = cp_firm;
            new_sup_cp = new_eq_cp(cp_firm);
            new_cp_adj(new_cust_cp,new_sup_cp) = 1;
        end
        
        %find avg shortest path distance
        new_cp_path_dist = avg_shortest_path_dist(new_cp_adj);
        %store it
        new_cp_conn_E(pw_index) = new_cp_path_dist;
        
        %find the index of that new eq network
        new_cp_eq_index = find(ismember(new_edge_list_cp,new_eq_cp,'rows'));
        
        %find the new agg output measure
        new_cp_output = new_agg_output_cp(new_cp_eq_index);
        
        %store the new output
        new_output_cp(pw_index) = new_cp_output;
        
        %calculate this relative output and store it
        cp_relative_output = new_cp_output/cp_output;
        cp_relative(pw_index) = cp_relative_output;
        
        %store the connectivity
        big_cp_connectivity_E(pw_index) = avg_dist_cp;
        
        %set cp flag to 1
        cp_flag_edge(pw_index) = 1;
            
        %%END SUCCESS CP EQ FIND CHECK    
        end
        
        
        %increment pw_index
        if pw_index < num_edges_deleted
            temp_pw_index = pw_index;
            pw_index = temp_pw_index + 1;
        end    
    end

    
    %for each edge in the planner solution, delete it, calculate new eq
    for pl_c = 1:num_firms
        %get the identity of the supplier for this edge in the original eq
        pl_edge_supplier = planner_network(pl_c);
        
        %record that edge's customer's possible suppliers and number of
        %customers in eq
        planner_poss_sup(pl_index) = num_pot_sups(pl_c);
        planner_num_cust(pl_index) = num_planner_cust(pl_c);
        
        %delete the agg outputs that correspond to that
        %particular edge
        
        %anywhere a pw_edge_supplier appears in the pw_c column, delete
        %that row in them all.
        %this is the list of rows to delete
        delete_index_list_pl = find(eq_edge_list_pl(:,pl_c)==pl_edge_supplier);
        
        %copy the output list
        new_ut = Ut;
        %cope the edge list
        new_edge_list_pl = eq_edge_list_pl;
        
        %for each row in the list to delete, delete that item from the new
        %utility list
        new_ut(delete_index_list_pl,:) = [];
        
        %for each row in the list of networks, deleted those 
        new_edge_list_pl(delete_index_list_pl,:) = [];
        
        %find the new max utility
        [new_planner_output, new_planner_index] = max(new_ut);
        
        %what is the associated network
        new_eq_pl = new_edge_list_pl(new_planner_index);
        
        %find the connectivity of the new network
        %first, turn it into an adjacency matrix
        new_planner_adj = zeros(num_firms);
        %for each firm, put a 1 in the place of it's supplier in the adjacency
        %matrix
        for frm_pl = 1:num_firms
            new_cust_pl = frm_pl;
            new_sup_pl = new_eq_pw(frm_pl);
            new_planner_adj(new_cust_pl,new_sup_pl) = 1;
        end
        %find avg shortest path distance
        new_pl_path_dist = avg_shortest_path_dist(new_planner_adj);
        %store it
        new_pl_conn_E(pl_index) = new_pl_path_dist;
        
        %how many customers does the pl_c have in the new network?
        %How many times does the new edge appear in the new eq network?
        new_num_cus_pl = sum(new_eq_pl == pl_c);
        %store it
        new_num_cust_pl(pl_index) = new_num_cus_pl;
        
        %store the new output
        new_output_pl(pl_index) = new_planner_output;
        
        %calculate this relative output and store it
        pl_relative_output = new_planner_output/planner_output;
        planner_relative(pl_index) = pl_relative_output;
        
        %store the connectivity
        big_planner_connectivity_E(pl_index) = avg_dist_pl;
        
        %set the planner flag to 1
        planner_flag_edge(pl_index) = 1;
        
        %increment the index
        if pl_index < num_edges_deleted
            temp_pl_index = pl_index;
            pl_index = temp_pl_index + 1;
        end
    end
    
    
    %calculate relative output
    rel_out = pwise_output/planner_output;
    %store it
    planner_relative_output(rep) = rel_out;
    
    %calculate the cp and planner relative output
    cp_rel_out = cp_output/planner_output;
    cp_planner_relative_output(rep) = cp_rel_out;
    
    %calculate the pwise and cp relative output
    pwise_cp_rel = pwise_output/cp_output;
    cp_pwise_relative_output(rep) = pwise_cp_rel;
    
    %create start and end indices for averages
    start = pl_index - num_firms+1;
    stop = pl_index;
    
    %calculate the average relative output in the planner experiment
    average_pl_output = mean(planner_relative(start:stop));
    avg_planner_relative(rep) = average_pl_output;
    
    %calculate the average relative output in the pwise experiment
    average_pw_output = mean(pwise_relative(start:stop));
    avg_pwise_relative(rep) = average_pw_output;
    
    %END FOR BOTH SUCCESSFUL
    end
end

%find the successful flag
success_flag = planner_flag_A.*pwise_flag_A;
%make the big successful flag
big_success_flag = planner_flag_edge.*pwise_flag_edge.*cp_flag_edge;


%how many successful experiments?
disp('How many successful originals?')
nnz(success_flag)
disp('How many successful edge deletion experiments?')
nnz(big_success_flag)

%Output Experiment One
disp('Experiment One: Planner vs. Pairwise vs. Coalition-Proof')

%planner_flag_A
%pwise_flag_A
%connectivity_A
%planner_relative_output

%get only the values where you successfully found a solution for the
%planner
good_connectivity = connectivity_A(success_flag == 1);
good_output = planner_relative_output(success_flag == 1);
good_cp_output = cp_planner_relative_output(success_flag == 1);
good_pwise_cp_output = cp_pwise_relative_output(success_flag == 1);

%big connectivity of A
big_A_conn = repelem(connectivity_A,5);
good_big_A_conn = big_A_conn(big_success_flag == 1);

%copy for trimming
perc_output = good_output;
trim_conn = good_connectivity;
trim_output = good_output;

%get good_output percentile
rel_p = prctile(perc_output, [1 99]);

%find the indices of output that lie outside the 1-99 percentile
p_del_1_below = find(perc_output <= rel_p(1));
p_del_1_above = find(perc_output >= rel_p(2));
%concatenate them
p_del_1 = cat(1, p_del_1_below, p_del_1_above);

%delete the values that correspond to those indices
trim_output(p_del_1) = [];
trim_conn(p_del_1) = [];

%regression (using trimmed data)
disp('Trimmed Pwise:')
[pl_pw_bt, pl_pw_intt] = regress(trim_output,trim_conn)
disp('Untrimmed Pwise:')
[pl_pw_b, pl_pw_int] = regress(good_output,good_connectivity)

%the cp planner one
disp('Coalition-Proof/Planner:')
[cp_rel_b, cp_rel_int] = regress(good_cp_output,good_connectivity)

%the pwise vs cp one
disp('Pairwise/Coalition-Proof:')
[pw_cp_rel_b, pw_cp_rel_int] = regress(good_pwise_cp_output, good_connectivity)

%Output Experiment Two
disp('Experiment Two: Planner Edge Deletion')
disp('Part I: Connectivity vs Average Relative Output')
%planner_flag_A
%connectivity_A
%planner_connectivity_E
%avg_planner_relative
%correlation_avg_planner_deleted = corrcoef(connectivity_A,avg_planner_relative)
%correlation_avg_planner_conn_deleted = corrcoef(planner_connectivity_E,avg_planner_relative)

%create good variables
good_pl_rel = planner_relative(big_success_flag == 1);
good_pl_conn_E = big_planner_connectivity_E(big_success_flag == 1);

%regression
[pl_c_b, pl_c_int] = regress(good_pl_rel, good_pl_conn_E)


disp('Part II: Centrality vs. Relative Output')
%planner_flag_edge
%planner_poss_sup
%planner_num_cust
%planner_relative

%create good ones
good_pl_poss_sup = planner_poss_sup(big_success_flag == 1);
good_pl_num_cust = planner_num_cust(big_success_flag == 1);

planner_X = [good_pl_poss_sup good_pl_num_cust];
[pl_e_b, pl_e_int] = regress(good_pl_rel, planner_X)

%Output Experiment Three
disp('Experiment Three: Pairwise Eq Edge Deletion')
disp('Part I: Connectivity vs Average Relative Output')
%pwise_flag_A
%connectivity_A
%pwise_connectivity_E
%avg_pwise_relative
%correlation_avg_pwise_deleted = corrcoef(connectivity_A,avg_pwise_relative)
%correlation_avg_pwise_conn_deleted = corrcoef(pwise_connectivity_E,avg_planner_relative)

%create good ones
good_pw_rel = pwise_relative(big_success_flag == 1);
good_pw_conn_E = big_pwise_connectivity_E(big_success_flag == 1);

%copy vars for trimming
perc_pw_rel = good_pw_rel;
trim_pw_rel = good_pw_rel;
trim_pw_conn = good_pw_conn_E;

%find 1, 99th percentile of pw rel
pw_p = prctile(perc_pw_rel, [1 99]);

%get the indices that fall outside the 1-99
p_del_3a_below = find(perc_pw_rel <= pw_p(1));
p_del_3a_above = find(perc_pw_rel >= pw_p(2));
%concatenate them
p_del_3a = cat(1, p_del_3a_below, p_del_3a_above);

%delete the associated values
trim_pw_rel(p_del_3a) = [];
trim_pw_conn(p_del_3a) = [];


%regression
disp('Trimmed:')
[pw_c_bt, pw_c_intt] = regress(trim_pw_rel, trim_pw_conn)
disp('Untrimmed:')
[pw_c_b, pw_c_int] = regress(good_pw_rel, good_pw_conn_E)



disp('Part II: Centrality vs. Relative Output')
%pwise_flag_edge
%pwise_poss_sup
%pwise_num_cust
%pwise_relative

%good ones
good_pw_poss_sup = pwise_poss_sup(big_success_flag == 1);
good_pw_num_cust = pwise_num_cust(big_success_flag == 1);
pwise_X = [good_pw_poss_sup good_pw_num_cust];

%copy for trim
trim_poss_sup = good_pw_poss_sup;
trim_num_cust = good_pw_num_cust;

%delete the associated values
trim_poss_sup(p_del_3a) = [];
trim_num_cust(p_del_3a) = [];


trim_pwise_X = [trim_poss_sup trim_num_cust];
disp('Trimmed:')
[pw_e_bt, pw_e_intt] = regress(trim_pw_rel, trim_pwise_X)
disp('Untrimmed:')
[pw_e_b, pw_e_int] = regress(good_pw_rel, pwise_X)

disp('Experiment Four: Coalition-Proof Eq Edge Deletion')
disp('Part I: Connectivity vs Average Relative Output')

%create good ones
good_cp_rel = cp_relative(big_success_flag == 1);
good_cp_conn_E = big_cp_connectivity_E(big_success_flag == 1);

%regression
[cp_rel_b, cp_rel_int] = regress(good_cp_rel, good_cp_conn_E)

disp('Part II: Centrality vs. Relative Output')

%create good ones
good_cp_poss_sup = cp_poss_sup(big_success_flag == 1);
good_cp_num_cust = cp_num_cust(big_success_flag == 1);

cp_X = [good_cp_poss_sup good_cp_num_cust];
[cp_c_b, cp_c_int] = regress(good_cp_rel, cp_X)


%Customer and Connectivity Output
disp('Customers and Connectivity')

%create good connectivities
good_big_conn_pl = big_planner_connectivity_E(big_success_flag == 1);
good_big_conn_pw = big_pwise_connectivity_E(big_success_flag == 1);
good_big_conn_cp = big_cp_connectivity_E(big_success_flag == 1);

%planner regression
[pl_cust_b, pl_cust_int] = regress(good_big_conn_pl, good_pl_num_cust)

%pwise regression
[pw_cust_b, pw_cust_int] = regress(good_big_conn_pw, good_pw_num_cust)

%cp regression
[cp_cust_b, cp_cust_int] = regress(good_big_conn_cp, good_cp_num_cust)

disp('Network Characteristics Regression:')
%create good new eq measures
good_cp_conn_E_new = new_cp_conn_E(big_success_flag == 1);
good_cp_num_cust_new = new_num_cust_cp(big_success_flag == 1);

disp('What makes it more likely for output to increase?')
%create indicator of output being over 1
ind_cp_rel = (good_cp_rel > 1);
ind_X = [good_big_A_conn good_cp_conn_E good_cp_conn_E_new good_cp_poss_sup good_cp_num_cust good_cp_num_cust_new];

[ind_b, ind_dev, ind_stats] = glmfit(ind_X,ind_cp_rel,'binomial','link','logit')
ind_stats.p

%calculate R^2
ind_residuals = ind_stats.resid;
ind_SSR = sum(ind_residuals.*ind_residuals);
ind_y_bar = mean(ind_cp_rel);
ind_y_diff = ind_cp_rel - ind_y_bar;
ind_SST = sum(ind_y_diff.*ind_y_diff);
ind_R_sq = 1 - (ind_SSR/ind_SST)

disp('Conditional on output increasing, what makes it increase by more?')
above_one_cp_out = good_cp_rel(ind_cp_rel == 1);
cond_good_big_A_conn = good_big_A_conn(ind_cp_rel == 1);
cond_good_cp_conn_E = good_cp_conn_E(ind_cp_rel == 1);
cond_good_cp_conn_E_new = good_cp_conn_E_new(ind_cp_rel == 1);
cond_good_cp_poss_sup = good_cp_poss_sup(ind_cp_rel == 1);
cond_good_cp_num_cust = good_cp_num_cust(ind_cp_rel == 1);
cond_good_cp_num_cust_new = good_cp_num_cust_new(ind_cp_rel == 1);

cond_X = [cond_good_big_A_conn cond_good_cp_conn_E cond_good_cp_conn_E_new cond_good_cp_poss_sup cond_good_cp_num_cust cond_good_cp_num_cust_new];


[cond_b, cond_int] = regress(above_one_cp_out, cond_X)


