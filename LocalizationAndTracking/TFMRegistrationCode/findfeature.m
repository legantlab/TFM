function featcoord = findfeature(beadpos,zmin_raw,zmax_raw)
scale_factor = [20,20,5];
beadpos_scaled = round(beadpos./scale_factor);
x_min = min(beadpos_scaled(:,1));
x_max = max(beadpos_scaled(:,1));
y_min = min(beadpos_scaled(:,2));
y_max = max(beadpos_scaled(:,2));
%z_min = min(beadpos_scaled(:,3));
z_min = floor(zmin_raw/scale_factor(3));
%z_max = max(beadpos_scaled(:,3));
z_max = ceil(zmax_raw/scale_factor(3));
beadpos_scaled = beadpos_scaled(beadpos_scaled(:,3)>=z_min,:);
beadpos_scaled = unique(beadpos_scaled,'rows');
featcoord = [];
featnum = 1;
for i = x_min:x_max
    for j = y_min:y_max
        for k = z_min:z_max
            if ~ isequal([i,j,k],beadpos_scaled(1,:))
                featcoord(featnum,:) = [i,j,k];
                featnum = featnum+1;
            else
                if length(beadpos_scaled(:,1))>1
                    beadpos_scaled = beadpos_scaled(2:end,:);
                else
                    return
                end
            end
        end
    end
end
