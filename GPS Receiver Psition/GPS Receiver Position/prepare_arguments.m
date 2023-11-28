function [satpos, obs_pseudo, weight] = prepare_arguments(broadcast_obs)
    satpos = zeros(3, size(broadcast_obs.data, 1));
    obs_pseudo = zeros(1, size(broadcast_obs.data, 1));
    weight = zeros(1, size(broadcast_obs.data, 1));
    for mm = 1: size(broadcast_obs.data, 1)
        for nn = 1:3
            satpos(nn, mm) = broadcast_obs.data(mm, nn);
        end
        obs_pseudo(mm) = broadcast_obs.data(mm, broadcast_obs.col.CorrP);
        weight_L1(mm) = broadcast_obs.data(mm, broadcast_obs.col.S1);
        weight_L2(mm) = broadcast_obs.data(mm, broadcast_obs.col.S2);
        weight(mm) = min(weight_L1(mm), weight_L2(mm));
    end
end