function ret = mirror_shards_distribute(n_run, n_shards)
% mirror_shards_distribute()
%
% Breaks a data array down into a number of shards, for use on single-node
% local worker pools, so we don't have to deal with Matlab DCE
% limitations on cores.
% Feed it a run number, and the number of shards to break the data over.
% The fundamentals of the simulation are all set here.

    q = 1;
    m = 1;
    nt = 10000;   % # timesteps
    dt = .01;     % step length
    qE = 0;
    qmt2 = q/m*dt/2;
    
    B0 = 50e-6; % Magnetic field base is 50 uT
    v0 = 0.00989179273; % likewise velocity base in PSL is equivalent to 25 eV
    r0 = 0.337212985; % based on Larmour radius w/ above, length base is ~0.337 m
    t0 = 7.14477319e-7; % based on B, Larmour period ~714 ns in s

    target_length = 5000;   % in km
    target_z = -target_length*1000/r0;  % negative because we're launching upwards
    long_enough = 1000000000;
    mirror_ratio = 5;
    saved_steps = 1000;

    % So Bsim=Breal/50uT, vsim=vreal/25 eV, and xsim=xreal/0.337m
    % So a 100x100x1000 simulation extent is a 33.7x33.7x337m volume
    % So dt ~71.4ns, and 1000 timesteps is 71us
    
    length_factor = target_z^2/(mirror_ratio-1); % assumes 'end point' is z=0
    
    x_range = 0;
    y_range = 0;
    z_range = target_z;
    v_range = [ 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, ...
        256, 289, 324, 361, 400, 441, 484, 529, 576, 625, 676, 729, ...
        784, 841, 900, 961, 1024, 1089, 1156, 1225 ]; % linear in v
    t_dphi = 3*pi/256; % delta for co-latitude
    t_domega = 0.001; % delta for solid angle in steradians
    %p_range = 0:pi/7:pi; %0:pi/15:pi/2;
    
    v_distrib = build_distrib(v0, x_range, y_range, z_range, v_range, t_dphi, t_domega);

    % Re-run particles that failed due to max timestep limit in Run 4
    load tzind.mat t_zind
    v_distrib = v_distrib(:,t_zind)

    parpool('torque_4nodes',n_shards)
    
    N_part = size(v_distrib,2);
    disp([ 'Distributing ' num2str(N_part) ' particles over ' num2str(n_shards) ' node shards...' ])
    
    v_sdivdist = distributed(v_distrib);
    
    disp('Start')
    tic

    % spmd (single program, multiple data) is a more generalized 
    % multithreaded methodology than parfor, and allows use of 
    % distributed/codistributed functionality to split up arrays
    spmd

        v_localdist = getLocalPart(v_sdivdist);
        N_dpart = size(v_localdist, 2);
        chunk_inds = globalIndices(v_sdivdist,2);
        
        disp( [ 'Shard ' num2str(labindex) ': ' num2str(N_dpart) ' particles (' num2str(chunk_inds(1)) ':' num2str(chunk_inds(end)) ').' ] )

        % have to call a function to use save inside an spmd.
        % because...raisins
        save_mah_data_plz(n_run, labindex, n_shards, N_dpart, chunk_inds, v_localdist);
        
    end
    
    toc
    disp('Done.')
    
    save(['mshards-r' num2str(n_run) '-master.mat'], ...
        'n_run', 'n_shards', 'N_part', 'v_distrib', ...
        'q', 'm', 'nt', 'dt', 'qE', 'qmt2', ...
        'B0', 'v0', 'r0', 't0', 'target_length', 'target_z', 'long_enough', ...
        'length_factor', 'mirror_ratio', 'saved_steps', 'v_range', 't_dphi');
    
    ret = 0;
    
end

function save_mah_data_plz(n_run, labindex, n_shards, N_shardpart, chunk_inds, v_sharddist)

    save(['mshard-r' num2str(n_run) '-' num2str(labindex) 'of' num2str(n_shards) '-input.mat'], ...
        'n_run', 'N_shardpart', 'chunk_inds', 'v_sharddist');

end

function d = build_distrib(v0, x_range, y_range, z_range, v_range, t_dphi, t_domega)
    % Build particle distribution

    % initial positions x y z
    % initial velocities v theta phi (magnitude, azimuth, co-latitude)
    %   mag 25:2000 eV, azi 0, el 0:pi/2
    % input as [ x y z v theta phi ] columns in v_distrib_raw

    t_phis = 0+t_dphi:t_dphi:pi/2-t_dphi; % range of phis, discard first (pole) and last (plane)
    
    angle_list = [ 0 0 ];
    for i=1:length(t_phis)
        t_phi = t_phis(i);
        angle_list = [ angle_list ; 0 t_phi ];
    end

    %angle_list = angle_list([ 1 2 3 4 19 20 21 22 38 39 40 41 ],:)
    
    % limit angles for tests
    %angle_list = angle_list(sin(4*angle_list(:,2)).^2 >= 0.995,:); % wedges in azimuthal angle
    
    v_distrib_raw = zeros(6,length(x_range)*length(y_range)*length(z_range)*length(v_range)*length(angle_list));
    vdr_ind = 1;
    for i=1:length(x_range)
        for j=1:length(y_range)
            for k=1:length(z_range)
                for l=1:length(angle_list)
                    for m=1:length(v_range)
                        v_distrib_raw(:,vdr_ind) = [ x_range(i) y_range(j) z_range(k) v_range(m) angle_list(l,1) angle_list(l,2) ];
                        vdr_ind = vdr_ind + 1;
                    end
                end
            end
        end
    end

    % Transform v_mag,theta,phi to v_x, v_y, v_z
    v_distrib = v_distrib_raw;
    t_v = sqrt(3.913903e-6*v_distrib_raw(4,:))/v0; % number is 2/(m_e*c^2) in eV^-1
    v_distrib(4,:) = t_v .* cos(v_distrib_raw(5,:)) .* sin(v_distrib_raw(6,:));
    v_distrib(5,:) = t_v .* sin(v_distrib_raw(5,:)) .* sin(v_distrib_raw(6,:));
    v_distrib(6,:) = t_v .*                            cos(v_distrib_raw(6,:));
    
    d = v_distrib;
end
