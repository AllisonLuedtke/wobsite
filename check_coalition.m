%this function checks to see if the given network has any profitable
%coaliton sized deviations

%it takes as input
%   -the network to check
%   -it's row number in the list of possible eq networks
%   -the list of possible eq networks
%   -the list of associated payoffs
%   -the list of suppliers for each firm
%   -the coalition size to check

function [dev] = check_coalition(Poss_Eq, Pot_index, Pot, Pay, Sup, coalition_size)

%get function variables
current = Poss_Eq;
current_index = Pot_index;
Pot_Max = Pot;
Payoffs = Pay;
Suppliers = Sup;
co_size = coalition_size;

%current payoffs
current_payoffs = Payoffs(current_index,:);
%initialize deviation as 0
dev = 0;

num_firms = length(current);

%make a list of the coalitions
coalition_list = combnk((1:num_firms),co_size);
%how many are there?
num_co = size(coalition_list,1);

%go through each coalition
co_count = 1;

%COALITIONS LOOP
while co_count <= num_co
    coalition = coalition_list(co_count,:);
    
    %create the list of other networks based on this coalition
    other_list = build_other_list(current, coalition, Suppliers);
    %how many are there?
    num_others = size(other_list,1);
    
    %peruse the other networks and check for profitable deviations
    other_count = 1;
    %OTHERS LOOP
    %if it's not empty ( = other_list(1,1) == 0)
    if other_list(1,1) ~= 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while other_count <= num_others
        %get the other network to compare to
        other_network = other_list(other_count,:);
        %the other suppliers to worry about are the suppliers of the
        %coalition firms in the other_network
        other_sups = other_network(coalition);
        %find the index of the other_network in the Pot Matrix
        other_index = find(ismember(Pot_Max,other_network,'rows'));
        %find the associated payoffs
        other_payoffs = Payoffs(other_index,:);
        
        %find the relevant payoffs
        %first, the payoffs of the coalition firms
        %before
        current_co_payoffs = current_payoffs(coalition);
        %after
        other_co_payoffs = other_payoffs(coalition);
        %then the payoffs of the new suppliers
        %before
        current_sup_payoffs = current_payoffs(other_sups);
        %after
        other_sup_payoffs = other_payoffs(other_sups);
        
        %check if BOTH coalition payoffs bigger after AND supplier payoffs
        %bigger after
        ones_vector = ones(1,co_size);
        %(other_co_payoffs > current_co_payoffs)
        %(other_sup_payoffs > current_sup_payoffs)
        if (isequal((other_co_payoffs > current_co_payoffs),ones_vector) && isequal((other_sup_payoffs > current_sup_payoffs),ones_vector))
        
            %you found a profitable deviation, set dev to 1
            dev = 1;
            %give up on this coalition, stop checking others
            other_count = num_others + 1;
        else
            %you didn't find a deviation
            %go to the next other network
            temp_other = other_count;
            other_count = temp_other + 1;
        end 
      
    %END OTHERS LOOP
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %END check if empty other list
    end
    %dev
    
    %if you found a deviation, give up on this coalition, return dev = 1
    if dev == 1
        %set coalition counter to end
        co_count = num_co + 1;
    else
        %you didn't find a deviation
        %you checked all the other networks for this particular coalition
        %go to the next coalition
        temp_co = co_count;
        co_count = temp_co + 1;
    end
    %if you get here and dev has never been switched to 1, you never found
    %a deviation and this coalition size had no profitable deviations  
    
%END COALITIONS LOOP    
end


end