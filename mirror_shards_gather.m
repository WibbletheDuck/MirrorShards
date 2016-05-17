function mirror_shards_gather(master_file)
% mirror_shards_gather()
%
% Compiles the output produced by per-node mirror shards that ran on input
% constructed by mirror_shards_distribute().

    % load the things we care about
    p_g = load(master_file, ...
        'n_run', 'n_shards', 'N_part', 'v_distrib', 'saved_steps');
    
    
    n_run = p_g.n_run; n_shards = p_g.n_shards; saved_steps = p_g.saved_steps;
    N_part = p_g.N_part; v_distrib = p_g.v_distrib;
    
    parpool('torque_4nodes',n_shards);
    
    disp([ 'Loading ' num2str(n_shards) ' shard outputs...' ])
    
    spmd
        
        [ v_sharddist, gm_result, gm_X, gm_V, gm_B ] = load_mah_data_plz(n_run, labindex, n_shards);
        
        t_codist_distrib = codistributor1d(2, codistributor1d.unsetPartition, [6, N_part]);
        t_codist_result = codistributor1d(3, codistributor1d.unsetPartition, [3, 7, N_part]);
        t_codist_saved = codistributor1d(3, codistributor1d.unsetPartition, [3, saved_steps, N_part]);
        
        % gather() to copy from GPU RAM to Main Memory
        r_divdist = codistributed.build(v_sharddist, t_codist_distrib, 'noCommunication');
        r_divres = codistributed.build(gm_result, t_codist_result, 'noCommunication');
        r_divsavX = codistributed.build(gm_X, t_codist_saved, 'noCommunication');
        r_divsavV = codistributed.build(gm_V, t_codist_saved, 'noCommunication');
        r_divsavB = codistributed.build(gm_B, t_codist_saved, 'noCommunication');

    end
    
    disp('Done, saving...')
    
    r_dist = gather(r_divdist);
    r_res = gather(r_divres);
    r_savX = gather(r_divsavX);
    r_savV = gather(r_divsavV);
    r_savB = gather(r_divsavB);
    
    if ~isequal(r_dist, v_distrib)
        disp('Rebuilt distribution does not equal OG distribution from master file!')
    end
    
    save([ 'mshards-r' num2str(n_run) '-final.mat' ], ...
        'n_run', 'n_shards', 'N_part', 'v_distrib', ...
        'r_dist', 'r_res', 'r_savX', 'r_savV', 'r_savB');
    
    disp('...great success?')
    
end

function [ l_dist, l_res, l_gm_X, l_gm_V, l_gm_B ] = load_mah_data_plz(i_n_run,labindex,n_shards)

    p_d = load([ 'mshard-r' num2str(i_n_run) '-' num2str(labindex) 'of' num2str(n_shards) '-output.mat' ], ...
        'r_shard_dist', 'r_shard_res', 'r_shard_X', 'r_shard_V', 'r_shard_B');
    
    l_dist = p_d.r_shard_dist;
    l_res  = p_d.r_shard_res;
    l_gm_X = p_d.r_shard_X;
    l_gm_V = p_d.r_shard_V;
    l_gm_B = p_d.r_shard_B;
    
end
