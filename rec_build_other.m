%This is the recursive function (called by build_other_list) that builds
%the list of other networks to check for a given current network and
%coalition.



function [others_mat] = rec_build_other(f, co, c_vec, b_vec, sup, o_mat)

%get function vars
firm_counter = f; %firm we're on
coalition = co; %coalition we're working with
current = c_vec;    %the current vector
build_vec = b_vec; %vector we're building to eventually store
Suppliers = sup;    %suppliers available to each firm
others_mat = o_mat; %matrix containing the list of others that we're building

%how many firms in the coalition?
co_size = length(coalition);

%RETURN CONDITION
%if we're on the last firm of the coalition, put it in the matrix listing
%them
if firm_counter == co_size
   %find the next available spot
   next_spot = find(others_mat(:,1)==0,1,'first'); 
   %store the existing build_vec in that spot
   others_mat(next_spot,:) = build_vec;
   
   
else
    %otherwise:
    co_firm = coalition(firm_counter+1);
    %find the suppliers of the next firm
    firm_sups = nonzeros(Suppliers(co_firm,:));
    %how many?
    num_firm_sups = length(firm_sups);
    %what's the current sup?
    current_sup = current(co_firm);
    for s = 1:num_firm_sups
        %new sup
        new_sup = firm_sups(s);
        %if it isn't the current sup, put it in the build vec for firm+1;
        if new_sup ~= current_sup
            
            %put it in as the supplier for firm + 1 and call the recursive function for the next firm
            build_vec(co_firm) = new_sup;
            others_mat = rec_build_other(firm_counter+1, coalition, current, build_vec, Suppliers, others_mat);
       end
        
    end
       
end









end