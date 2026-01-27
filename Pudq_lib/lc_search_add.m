function loop_closures = lc_search_add(loop_closures, i, j)
    
    pair_exists = any(loop_closures(1,:) == i & loop_closures(2,:) == j);

    if ~pair_exists
        loop_closures(:, end+1) = [i; j];
    end

end