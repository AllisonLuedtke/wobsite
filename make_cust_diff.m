function [diff_cust] = make_cust_diff(old_cust, new_cust)

%get function variables
old_cust_list = old_cust;
new_cust_list = new_cust;

%size?
num_spots = size(old_cust_list,1)*size(old_cust_list,2);
num_rows = size(old_cust_list,1);

%create blank diff cust list
diff_cust = zeros(num_spots,1);
count = 1;


for row = 1:num_rows
    %if there are customers
    if old_cust_list(row,1) ~= 0
        %for each nonzero element
        num_cust = nnz(old_cust_list(row,:));
        for c = 1:num_cust
            %if not infinity
            if (old_cust_list(row,c) ~= Inf) && (new_cust_list(row,c) ~= Inf)
                %find difference
                diff = old_cust_list(row,c) - new_cust_list(row,c);
                %store difference
                diff_cust(count) = diff;
                %increment count
                temp_count = count;
                count = temp_count + 1;
            end     
                  
        end
       
    end
        
end

%trim the diff_cust list
end_index = nnz(diff_cust);
diff_cust(end_index+1:end) = [];
       
end




