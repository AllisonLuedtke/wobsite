%Practicing the recursive function

Pot_Adj = [0 1 1; 1 0 1; 1 1 0];

num_firms = 3;

sup_matrix = [2 3; 1 3; 1 2];

init_vec = [0 0 0];

storing_matrix = zeros(8,3);


%start at firm 0
rec_matrix(0,init_vec,sup_matrix,storing_matrix)



%You'll eventually need to build the matrices from the lists of suppliers