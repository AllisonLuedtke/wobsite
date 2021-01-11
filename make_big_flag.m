%this function takes the flag that indicates the succesful reps (both
%planner and pwise were successful) and makes a big flag for each firm
%level (edge deletion) outcome
function [big_flag] = make_big_flag(flag, num_f)

%get function vars
suc_flag = flag;
num_firms = num_f;

old_flag_length = length(suc_flag);

%initialize big flag
big_flag = ones(old_flag_length*num_firms,1);

%find the zeros in the original flag
zero_list = find(suc_flag == 0);

%create a list of num_firms zeros to put in wherever you need them
zero_column = zeros(num_firms,1);

%for each zero, put 6 zeros in the big flag in the right place
for z = 1:length(zero_list)
    %get that index
    z_index = zero_list(z);
    %start index
    start_index = num_firms*(z_index - 1) + 1;
    %stop index
    stop_index = start_index + (num_firms - 1);
    %replace start-stop of the big flag with the zeros
    big_flag(start_index:stop_index) = zero_column;
    
end





end