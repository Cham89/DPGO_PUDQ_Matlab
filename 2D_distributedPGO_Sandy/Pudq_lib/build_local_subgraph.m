function local_graph = build_local_subgraph(multi_graph, r)
    
    local_graph.robot_id            = r;
    local_graph.anchor_first        = (r == 1); 
    local_graph.vertices_true       = multi_graph(r).vertices_true;
    local_graph.vertices_pudq_true  = multi_graph(r).vertices_pudq_true;
    local_graph.vertices            = multi_graph(r).vertices;
    local_graph.vertices_pudq       = multi_graph(r).vertices_pudq;
    local_graph.vertices_interval = multi_graph(r).vertices_interval;

    local_graph.lm_vertices         = multi_graph(r).lm_vertices;
    local_graph.lm_vertices_pudq    = multi_graph(r).lm_vertices_pudq;
    local_graph.lm_vertices_true    = multi_graph(r).lm_vertices_true;
    local_graph.lm_vertices_pudq_true = multi_graph(r).lm_vertices_pudq_true;
    local_graph.lm_foreign_info      = multi_graph(r).lm_foreign_info;

    local_graph.intra_edges         = multi_graph(r).intra_edges;        
    local_graph.intra_dp            = multi_graph(r).intra_dp;
    local_graph.intra_dp_true       = multi_graph(r).intra_dp_true;
    local_graph.intra_dp_pudq       = multi_graph(r).intra_dp_pudq;
    local_graph.intra_dp_pudq_true  = multi_graph(r).intra_dp_pudq_true;
    local_graph.intra_t             = multi_graph(r).intra_t;
    local_graph.intra_R             = multi_graph(r).intra_R;
    local_graph.intra_info          = multi_graph(r).intra_info;
    local_graph.intra_info_pudq     = multi_graph(r).intra_info_pudq;
    local_graph.intra_info_se2      = multi_graph(r).intra_info_se2;
    local_graph.intra_tau           = multi_graph(r).intra_tau;
    local_graph.intra_kappa         = multi_graph(r).intra_kappa;

    local_graph.inter_edges         = multi_graph(r).inter_edges;         
    local_graph.inter_dp            = multi_graph(r).inter_dp;
    local_graph.inter_dp_true       = multi_graph(r).inter_dp_true;
    local_graph.inter_dp_pudq       = multi_graph(r).inter_dp_pudq;
    local_graph.inter_dp_pudq_true  = multi_graph(r).inter_dp_pudq_true;
    local_graph.inter_t             = multi_graph(r).inter_t;
    local_graph.inter_R             = multi_graph(r).inter_R;
    local_graph.inter_info          = multi_graph(r).inter_info;
    local_graph.inter_info_pudq     = multi_graph(r).inter_info_pudq;
    local_graph.inter_info_se2      = multi_graph(r).inter_info_se2;
    local_graph.inter_tau           = multi_graph(r).inter_tau;
    local_graph.inter_kappa         = multi_graph(r).inter_kappa;

    local_graph.Omega_ij          = multi_graph(r).Omega_ij;

end
