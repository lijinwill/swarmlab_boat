function [vel_command, collisions] = compute_vel_vasarhelyi(self, p_swarm, r_agent, dt)
    % VASARHELYI SWARM ALGORITHM
    % This is an implementation of the Vasarhelyi algorithm. It allows the
    % navigation of a swarm of agents in presence of obstacles and walls.
    %
    % Ref:      Vasarhelyi, Science Robotics, 2018
    %       Vásárhelyi G, Virágh C, Somorjai G, et al. Optimized flocking of autonomous drones
    %       in confined environments[J]. Science Robotics, 2018, 3(20): eaat3536.
    % 
    % Modif:    a cohesion term has been added to make the agents get
    %           closer when they are farther than r0_rep.
    %
    % Inputs:
    %   p_swarm: swarm parameters
    %   r_agent: safety radius of agents
    %   dt: time step
    %
    % Outputs:
    %   vel_command: commanded velocities for every agent
    %   collisions: [nb_agent_collisions nb_obs_collisions min_dist_obs]
    %
    
    
    %% Initialize variables

    pos = self.get_pos_ned();
    vel = self.get_vel_ned();

    % Initialize variables
    nb_agents = self.nb_agents;
    M = zeros(nb_agents, nb_agents);    % Neighborhood matrix
    dist_mat = zeros(nb_agents, nb_agents);    % Distance matrix
    vel_rep = zeros(3, nb_agents);      % Repulsion velocity
    vel_fric = zeros(3, nb_agents);     % Velocity matching velocity
    vel_wall = zeros(3, nb_agents);     % Arena repulsion velocity
    vel_obs = zeros(3, nb_agents);      % Obstacle repulsion velocity
    vel_command = zeros(3, nb_agents);  % Total commanded velocity
    
    nb_agent_collisions = 0; % Nb of collisions among agents
    nb_obs_collisions = 0; % Nb of collisions against obstacles
    min_dist_obs = 20; % 猜测是离障碍物的最小距离

    
    %% Compute velocity commands
    
    for agent = 1:nb_agents
        
        
        %% Find neighbors
        
        % Compute agent-agent distance matrix
        p_rel = pos - pos(:, agent);% 其他agent坐标减去当前agent坐标
        dist = sqrt(sum((p_rel.^2), 1));
        dist_mat(agent, :) = dist;

        % Define neighbours list
        neig_list = (1:nb_agents)';
        neig_list = neig_list(dist ~= 0);

        % Count collisions
        nb_agent_collisions = nb_agent_collisions + sum(dist < 2 * r_agent) - 1;

        % Initialize number of neighbours
        nb_neig = nb_agents - 1;

        % Constraint on neighborhood given by the euclidean distance
        % 这里计算邻居数量nb_neig，其中p_swarm.r来自param_swarm，表示最大影响距离，小于此距离的都是邻居
        if isfield(p_swarm, 'r')
            neig_list = neig_list(dist(neig_list) < p_swarm.r);
            nb_neig = length(neig_list);
        end

        % Constraint on neighborhood given by the topological distance
        if isfield(p_swarm, 'max_neig')
            if nb_neig > p_swarm.max_neig
                [~, idx] = sort(dist(neig_list));
                neig_list = neig_list(idx(1:p_swarm.max_neig));
                nb_neig = p_swarm.max_neig;
            end
        end

        % Adjacency matrix (asymmetric in case of limited fov)
        % 矩阵形式记录邻居，下面再也没用过
        M(agent, neig_list) = 1;

        
        %% Compute different contributions
        
        if nb_neig ~= 0
            v_rel = vel - vel(:, agent);
            v_rel_norm = sqrt(sum((v_rel.^2), 1));%类似于距离dist，只不过是用速度矢量算出来的

            % Compute vel and pos unit vector between two agents
            % 这个p_rel_u就是论文里面的ri-rj/rij，表明一种方向性
            p_rel_u = -p_rel ./ dist;
            v_rel_u = -v_rel ./ v_rel_norm;

            for agent2 = neig_list'
                
                % Repulsion and attraction
                if dist(agent2) < p_swarm.r0_rep  % repulsion
                    %下面这行代码来自论文第7页和第8页的公式2和3
                    vel_rep(:, agent) = vel_rep(:, agent) + ...
                        p_swarm.p_rep * (p_swarm.r0_rep - dist(agent2)) * p_rel_u(:, agent2);
                else  % attraction
                    vel_rep(:, agent) = vel_rep(:, agent) + ...
                        p_swarm.p_rep * (dist(agent2) - p_swarm.r0_rep) *- p_rel_u(:, agent2);
                end

                % Velocity alignement
                % 公式4和公式5
                v_fric_max = get_v_max(p_swarm.v_fric, dist(agent2) - p_swarm.r0_fric, p_swarm.a_fric, p_swarm.p_fric);

                % 公式6
                if v_rel_norm(agent2) > v_fric_max
                    vel_fric(:, agent) = vel_fric(:, agent) + ...
                        p_swarm.C_fric * (v_rel_norm(agent2) - v_fric_max) * v_rel_u(:, agent2);
                end
            end
        end
        
        
        %% Wall and obstacle avoidance
        % 下面三个模块的并列的，第一个是关于边界墙的，第三个是有圆柱形障碍物的
        % Add arena repulsion effect
        if (p_swarm.is_active_arena == true)
            unit = eye(3);
            %On each axis we have the two repulsions
            for axis = 1:3
                %On each axis there is two forces (each side of the arena)
                for dir = 1:2
                    dist_ab = abs(pos(axis, agent) - p_swarm.x_arena(axis, dir));

                    %Compute velocity of wall shill agent toward center of the arena
                    v_wall_virtual = unit(:, axis) .* p_swarm.v_shill;

                    if dir == 2
                        v_wall_virtual = -v_wall_virtual;
                    end

                    %Compute relative velocity (Wall - Agent)
                    vel_ab = sqrt(sum((vel(:, agent) - v_wall_virtual).^2));

                    v_wall_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                    if vel_ab > v_wall_max
                        vel_wall(:, agent) = vel_wall(:, agent) + ...
                            (vel_ab - v_wall_max) * (v_wall_virtual - vel(:, agent)) ./ vel_ab;
                    end
                end
            end
        end

        % Compute spheric effect
        if (p_swarm.is_active_spheres == true)

            for obs = 1:p_swarm.n_spheres
                % Get obstacle center and radius
                c_obs = p_swarm.spheres(1:3, obs);
                r_obs = p_swarm.spheres(4, obs);

                % Compute distance agent(a)-obstacle(b)
                dist_ab = sqrt(sum((pos(:, agent) - c_obs).^2)) - r_obs;
                nb_obs_collisions = nb_obs_collisions + sum(dist_ab < r_agent);

                % Set the virtual speed of the obstacle direction out of
                % the obstacle
                v_obs_virtual = (pos(:, agent) - c_obs) / (dist_ab + r_obs) * p_swarm.v_shill;

                % Compute relative velocity agent-obstacle
                vel_ab = sqrt(sum((vel(:, agent) - v_obs_virtual).^2));

                if dist_ab < min_dist_obs
                    min_dist_obs = dist_ab;
                end

                v_obs_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                if vel_ab > v_obs_max
                    vel_obs(:, agent) = vel_obs(:, agent) + (vel_ab - v_obs_max) * (v_obs_virtual - vel(:, agent)) ./ vel_ab;
                end
            end
        end

        % Compute cylindric effect
        if (p_swarm.is_active_cyl == true)

            for obs = 1:p_swarm.n_cyl
                % Get obstacle center and radius
                c_obs = p_swarm.cylinders(1:2, obs);
                r_obs = p_swarm.cylinders(3, obs);

                % Compute distance agent(a)-obstacle(b)
                dist_ab = sqrt(sum((pos(1:2, agent) - c_obs).^2)) - r_obs;
                nb_obs_collisions = nb_obs_collisions + sum(dist_ab < r_agent);

                % Set the virtual speed of the obstacle direction out of
                % the obstacle
                v_obs_virtual = (pos(1:2, agent) - c_obs) / (dist_ab + r_obs) * p_swarm.v_shill;

                % Compute relative velocity agent-obstacle
                vel_ab = sqrt(sum((vel(1:2, agent) - v_obs_virtual).^2));

                if dist_ab < min_dist_obs
                    min_dist_obs = dist_ab;
                end
                
                % 第8页公式8
                v_obs_max = get_v_max(0, dist_ab - p_swarm.r0_shill, p_swarm.a_shill, p_swarm.p_shill);

                % 第8页公式9
                if vel_ab > v_obs_max
                    vel_obs(1:2, agent) = vel_obs(1:2, agent) + ...
                        (vel_ab - v_obs_max) * (v_obs_virtual - vel(1:2, agent)) ./ vel_ab;
                end

            end

        end
        
        %% Sum agent-agent and obstacle contributions

        % 公式10
        vel_command(:, agent) = vel_rep(:, agent) + vel_fric(:, agent) + vel_obs(:, agent) + vel_wall(:, agent);

        % Add self propulsion OR migration term
        % 下面代码的意思是公式10最左边那个项，即给予的集群基础速度的形式与大小
        % 可以是向着一个方向迁徙，也可以是朝着一个给定的目标
        v_norm = sqrt(sum((vel(:, agent).^2), 1));

        if p_swarm.is_active_migration == true% migration
            vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * p_swarm.u_ref;
        elseif p_swarm.is_active_goal == true
            x_goal_rel = p_swarm.x_goal - pos(:, agent);
            u_goal = x_goal_rel / norm(x_goal_rel);
            vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * u_goal;
        else
            % self-propulsion，即每个个体没有方向漫无目的的运动
            if v_norm > 0
                vel_command(:, agent) = vel_command(:, agent) + p_swarm.v_ref * vel(:, agent) / v_norm;
            end
        end
    end
    
    
    %% Compute collisions and bound velocities and accelerations

    % Total number of collisions per time step
    nb_agent_collisions = nb_agent_collisions / 2; % reciprocal
    collisions = [nb_agent_collisions nb_obs_collisions min_dist_obs];

    % Add random effect on velocities
    if isfield(p_swarm, 'c_r')
        vel_command = vel_command + p_swarm.c_r * randn(3, nb_agents);
    end

    % Bound velocities and acceleration
    % 这个部分的作用是限制最大速度，论文公式11
    if ~isempty(p_swarm.max_v)
        vel_cmd_norm = sqrt(sum((vel_command.^2), 1));
        %v_norm = sqrt(sum((vel.^2), 1));
        
        idx_to_bound = (vel_cmd_norm > p_swarm.max_v);%这个表明此轮速度中需要被限制的agent的idx
        % 如果存在这种个体，则将速度命令修改为最大速度*该速度原先的方向
        if sum(idx_to_bound) > 0
            vel_command(:, idx_to_bound) = p_swarm.max_v * ...
                vel_command(:, idx_to_bound) ./ repmat(vel_cmd_norm(idx_to_bound), 3, 1);
        end
    end
    % 类似上面，但是是限制最大加速度
    if ~isempty(p_swarm.max_a)
        accel_cmd = (vel_command-vel)./dt;
        accel_cmd_norm = sqrt(sum(accel_cmd.^2, 1));
        idx_to_bound = ( accel_cmd_norm > p_swarm.max_a | accel_cmd_norm < - p_swarm.max_a);
        if sum(idx_to_bound) > 0
            vel_command(:, idx_to_bound) = vel(:, idx_to_bound) + ...
                dt*p_swarm.max_a * accel_cmd(:, idx_to_bound) ./ ...
                repmat(accel_cmd_norm(idx_to_bound), 3, 1);
        end
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate V fric max
% 这function是论文第八页公式4和公式5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ v_fricmax ] = get_v_max(v_fric, r, a, p)

    if r < 0
        v_fricmax = 0;
    elseif r * p > 0 && r * p < a / p
        v_fricmax = r * p;
    else
        v_fricmax = sqrt(2 * a * r - a^2 / p^2);
    end

    if v_fricmax < v_fric
        v_fricmax = v_fric;
    end
end
