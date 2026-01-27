function [Lm_idx, Exists] = findOrAddLandmark(multi_graph, r_j, v_j_local)
    Exists = false;
    for lm = 1:length(multi_graph.lm_foreign_info)
        info = multi_graph.lm_foreign_info{lm};
        if info.robot == r_j && info.vertex == v_j_local
            Lm_idx = lm;
            Exists = true;
            return;
        end
    end
    Lm_idx = length(multi_graph.lm_foreign_info) + 1;
end