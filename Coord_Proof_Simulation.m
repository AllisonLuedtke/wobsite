

%Number of firms
num_firms = 5;

%Create an Adjacency Matrix and a Productivity Matrix
[A,Z] = Create_A_and_Z(num_firms);

%Find the Coordination Proof Equilibrium Networks
[Output, Eq_Network, residuals, x_star, Payoffs, Eq_Mat_List, Other_Agg_Outputs, Eq_Efficiency, Efficiencies ] = Coalition_Proof_Simulation(A,Z);