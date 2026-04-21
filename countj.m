function csj = countj(j, ss)

csj =0;
for i=1:length(ss)
    if ss(i) == j
        csj = csj+1;
    end
end

end