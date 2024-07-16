%real_full_fast_beam.m
for level = 1:num-1
    recov = check_same_buffer(sets_buffer(level,:),sets_buffer_fb(level,:),bw);
    if recov>0
        disp(level)
        disp(recov)
    end
end
for level = 1:num-1
    recov = max(abs(sort(log_like_buffer(level,:))-sort(log_like_buffer_fb(level,:))));

    if recov>0.0001
        disp(level)
        disp(recov)
    end

end

