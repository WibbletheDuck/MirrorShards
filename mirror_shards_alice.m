function ret = mirror_shards_alice(n_shard, n_cores, master_file)
% mirror_shards_alice()
%
% Does the actual simulation work in the mirror_shards_* system.

    p_g = load(master_file,'n_run','n_shards', ...
        'dt','qE','qmt2','r0', 'long_enough', ...
        'target_length','mirror_ratio', 'saved_steps');
    
    target_z = -p_g.target_length*1000/p_g.r0;  % negative because we're launching upwards
    length_factor = target_z^2/(p_g.mirror_ratio-1); % assumes 'end point' is z=0

    % load values from p_g
    n_run = p_g.n_run; n_shards = p_g.n_shards; dt = p_g.dt; qE = p_g.qE;
    qmt2 = p_g.qmt2; long_enough = p_g.long_enough; saved_steps = p_g.saved_steps;

    % n_run, N_shardpart, chunk_inds, v_sharddist;
    p_s = load(['mshard-r' num2str(n_run) '-' num2str(n_shard) 'of' num2str(n_shards) '-input.mat'], ...
        'N_shardpart','v_sharddist');

    N_shardpart = p_s.N_shardpart;
    v_sharddist = p_s.v_sharddist;
    
    parpool('local', n_cores);

    % This is local to the shard now, but we'll just redefine below for the
    % part of the distribution that's local to each worker.
    v_sdivdist = distributed(v_sharddist);
    
    disp([ 'Simulating ' num2str(N_shardpart) ' particles over INFINITE timesteps...' ])
    tic
    disp('Start')

    % spmd (single program, multiple data) is a more generalized 
    % multithreaded methodology than parfor, and allows use of 
    % distributed/codistributed functionality to split up arrays
    spmd

        v_localdist = getLocalPart(v_sdivdist);
        N_dpart = size(v_localdist, 2);
        chunk_inds = globalIndices(v_sdivdist,2);
        
        %d = gpuDevice();
        disp( [ 'Running ' num2str(N_dpart) ' particles (' num2str(chunk_inds(1)) ':' num2str(chunk_inds(end)) ') in Lab ' num2str(labindex) '.' ] )

        % Pre-allocate result arrays on GPU
        gm_X = zeros([ 3, saved_steps, N_dpart ], 'double'); % x,y,z
        gm_V = zeros([ 3, saved_steps, N_dpart ], 'double'); % vx,vy,vz
        gm_Bv = zeros([ 3, saved_steps, N_dpart ], 'double'); % Bx,By,Bz

        % redundant array of results, seven 3-vectors containing
        % 1,2 position and velocity at target-z (z_t)
        % 3, # of timestep before z_t, after, and actual calculated crossing time
        % 4,5 position and velocity of pre-z_t timestep
        % 6,7 position and velocity of post-z_t timestep
        gm_result = zeros([ 3, 7, N_dpart ], 'double');

        %d = gpuDevice();
        %t_tmem = d.TotalMemory;
        %t_umem = t_tmem-d.AvailableMemory;
        %disp( [ 'Memory Used: ' num2str(t_umem/1e9) '/' num2str(t_tmem/1e9) 'GB (' num2str(t_umem/t_tmem*100) '%)' ]);
        
        active_indices = 1:N_dpart;

        % Get B at initial positions
        % permute() lets us slot a (3,N) data peg into a (3,M,N) hole
        gm_X(:,end-1,:) = permute(v_localdist(1:3,:),[1 3 2]);
        gm_V(:,end-1,:) = permute(v_localdist(4:6,:),[1 3 2]);
        
        % Recall all arrays are (dimension, timestep, particles)
        gm_B_x = squeeze(-gm_X(1,end-1,active_indices).*gm_X(3,end-1,active_indices)/length_factor);
        gm_B_y = squeeze(-gm_X(2,end-1,active_indices).*gm_X(3,end-1,active_indices)/length_factor);
        gm_B_z = squeeze(1+gm_X(3,end-1,active_indices).^2/length_factor);
        
        % Calculate 2nd position with Boris Mover 
        gm_v_mh = squeeze(gm_V(:,end-1,:));
        gm_v_minus = gm_v_mh + qE;

        gm_B = [ gm_B_x gm_B_y gm_B_z ].'; gm_Bv(:,end,:) = gm_B;
        gm_t_vec = qmt2*gm_B;
        gm_s_vec = 2*gm_t_vec./(1+gm_t_vec.^2);
        gm_v_prime = gm_v_minus + cross(gm_v_minus,gm_t_vec,1);
        gm_v_plus = gm_v_minus + cross(gm_v_prime,gm_s_vec,1);

        gm_V(:,end,:) = 0.5 .* (gm_v_mh + gm_v_plus + qE);
        gm_X(:,end,:) = gm_X(:,end-1,:) + gm_V(:,end-1,:) .* dt;

        tstep = 1;
        % Loop until all particles are done, or we've 
        % done an absurd number of timesteps.
        while ~isempty(active_indices)
            tstep = tstep + 1;
            
            % shift saved-data matrices down one row
            gm_X(:,1:end-1,active_indices) = gm_X(:,2:end,active_indices);
            gm_V(:,1:end-1,active_indices) = gm_V(:,2:end,active_indices);
            gm_Bv(:,1:end-1,active_indices) = gm_Bv(:,2:end,active_indices);
            
            if labindex == 1 && mod(tstep,10000) == 0
                display(['Step ' num2str(tstep) ', ' num2str(length(active_indices)) ... 
                    ' particles active, min/max z = ' num2str(min(gm_X(3,end,active_indices))) '/' num2str(max(gm_X(3,end,active_indices))) '.'])
            end
            
            % Recall all arrays are (dimension, timestep, particles)
            gm_B_x = squeeze(-gm_X(1,end-1,active_indices).*gm_X(3,end-1,active_indices)/length_factor);
            gm_B_y = squeeze(-gm_X(2,end-1,active_indices).*gm_X(3,end-1,active_indices)/length_factor);
            gm_B_z = squeeze(1+gm_X(3,end-1,active_indices).^2/length_factor);

            gm_v_minus = squeeze(gm_V(:,end-1,active_indices) + qE);    % half-step due to E-field
            
            gm_B = [ gm_B_x gm_B_y gm_B_z ].'; gm_Bv(:,end,active_indices) = gm_B;
            gm_t_vec = qmt2*gm_B;
            gm_s_vec = 2*gm_t_vec./(1+gm_t_vec.^2);
            gm_v_prime = gm_v_minus + cross(gm_v_minus,gm_t_vec);   % these calculate the
            gm_v_plus = gm_v_minus + cross(gm_v_prime,gm_s_vec);    % B-field effects
            
            gm_V(:,end,active_indices) = gm_v_plus + qE;  % second half-step from E-field
            gm_X(:,end,active_indices) = gm_X(:,end-1,active_indices) + gm_V(:,end-1,active_indices) .* dt;
            
            % check if next z-pos passes the target plane z=0
            strike_indices = active_indices(gm_X(3,end,active_indices) >= 0); 
            if ~isempty(strike_indices)
                %display([ 'Timestep ' num2str(tstep) ': ' num2str(length(strike_indices)) ' strikes.' ]);
                % interpolate absolute strike XVT
                % NB: assumes 'target z' is z=0 plane
                t_nStrikes = length(strike_indices);
                t_V0 = squeeze(gm_V(:,end-1,strike_indices)); t_V1 = squeeze(gm_V(:,end,strike_indices)); % initial and final velocities
                t_X0 = squeeze(gm_X(:,end-1,strike_indices)); t_X1 = squeeze(gm_X(:,end,strike_indices)); % initial and final positions
                
                t_a01 = ( t_V1-t_V0 )/dt; % acceleration from x_pre to x_post
                
                t_t0t = ( -t_V0(3,:) + sqrt(t_V0(3,:).^2 - 2*t_a01(3,:).*t_X0(3,:)) )./t_a01(3,:); % time to z=0
                %display(['foo: ' num2str( size(t_V0) ) ' - ' num2str( size(t_a01) ) ' - ' num2str( size(t_t0t) ) '.'])
                t_Vt = t_V0 + bsxfun(@times,t_a01,t_t0t); % velocity at z=0
                t_Xt = t_X0 + bsxfun(@times,t_V0,t_t0t) + 0.5*bsxfun(@times,t_a01,t_t0t.^2); % complete position at z=0
                
                % target x, target v, times, x0, v0, x1, v1
                gm_result(:,1,strike_indices) = squeeze(t_Xt);
                gm_result(:,2,strike_indices) = squeeze(t_Vt);
                
                %display(['foo: ' num2str( size((tstep+t_t0t)*dt) ) ' - ' num2str( size(repmat(tstep,1,t_nStrikes)) ) ' - ' num2str( size(t_t0t) ) '.'])
                t_t = [ (tstep+t_t0t)*dt ; repmat(tstep,1,t_nStrikes) ; t_t0t ];
                %display(['bar: ' num2str( size(t_t) ) ' - ' num2str( size(gm_result) ) '.']);
                
                gm_result(:,3,strike_indices) = t_t;

                gm_result(:,4,strike_indices) = squeeze(t_X0);
                gm_result(:,5,strike_indices) = squeeze(t_V0);
                
                gm_result(:,6,strike_indices) = squeeze(t_X1);
                gm_result(:,7,strike_indices) = squeeze(t_V1);
                
                active_indices = active_indices(~ismember(active_indices,strike_indices));
            end

        end

	display(['Final timesteps: ' num2str(tstep) '.'])

        t_codist_result = codistributor1d(3, codistributor1d.unsetPartition, [3, 7, N_shardpart]);
        t_codist_saved = codistributor1d(3, codistributor1d.unsetPartition, [3, saved_steps, N_shardpart]);
        
        % build codist arrays
        r_divres = codistributed.build(gm_result, t_codist_result, 'noCommunication');
        r_divsavX = codistributed.build(gm_X, t_codist_saved, 'noCommunication');
        r_divsavV = codistributed.build(gm_V, t_codist_saved, 'noCommunication');        
        r_divsavB = codistributed.build(gm_Bv, t_codist_saved, 'noCommunication');
        
    end % spmd
    
    % gather() to recombine distributed arrays
    r_shard_res = gather(r_divres);
    r_shard_X = gather(r_divsavX);
    r_shard_V = gather(r_divsavV);
    r_shard_B = gather(r_divsavB);
    r_shard_dist = gather(v_sdivdist);

    save(['mshard-r' num2str(n_run) '-' num2str(n_shard) 'of' num2str(n_shards) '-output.mat'], ...
        'N_shardpart', 'r_shard_dist', ...
        'r_shard_res', 'r_shard_X', 'r_shard_V', 'r_shard_B');
    
    disp('End')
    toc

    ret = 0;
    
end
