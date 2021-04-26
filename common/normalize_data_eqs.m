function data = normalize_data_eqs(data,sz)

count = 0;
for iii = 1:length(sz)
    data((1:sz(iii))+count) = data((1:sz(iii))+count)/max(abs(data((1:sz(iii))+count) ));
    count = count + sz(iii);
end
