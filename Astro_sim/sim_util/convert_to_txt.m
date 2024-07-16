function [] =  convert_to_txt(filename)
C = who('-file',strcat(filename,'.mat'));

for i = 1:length(C)
    disp(C{i})
    disp(i)
    datastruct = load(filename, C{i});
    if strcmp(class(datastruct.(C{i})),'double')
        writematrix(datastruct.(C{i}),strcat(filename,C{i},'.txt'))

    elseif strcmp(class(datastruct.(C{i})),'cell') 
        sz = size(datastruct.(C{i}));
        if sz(1)==1

            for jj = 1:length(datastruct.(C{i}))
                fid = fopen(strcat(filename,C{i},'.txt'), 'at'); 
                fprintf(fid, num2str(jj));
                fclose(fid);
                writecell(datastruct.(C{i})(jj),strcat(filename,C{i},'.txt'),'WriteMode','append')
            
    
            end

        else
            for jj = 1:length(datastruct.(C{i}))
                for kk = 1:length(datastruct.(C{i}){jj})
                    fid = fopen(strcat(filename,C{i},num2str(jj),'.txt'), 'at'); 
                    fprintf(fid, num2str(kk));
                    fclose(fid);
                    writecell(datastruct.(C{i}){jj}(kk),strcat(filename,C{i},num2str(jj),'.txt'),'WriteMode','append')
              
                end
            end
        end
    else
        continue;
    end
end