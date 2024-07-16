%   real_full_fast_beam.m
for ii = 73:num
    disp(num2str(ii))
    state_buffer = sets_all_buffer(ii,:);
    if check_duplicate_buffer(state_buffer,100)
        disp(strcat("duplicate ", num2str(ii)))
    end
end