function [] = plot_err(rec_xyz_actual, epochs, varargin)

nepoch = size(epochs,1);

    for i = 1 : (nargin-2)
        userPos = varargin{i};
        a = inputname(i+2);
        a = a(4:15);
        method_name(i,:) = a;
        for j = 1 : size(userPos,1)
            err_x(j) = abs(userPos(j,1) - rec_xyz_actual(1));
            err_y(j) = abs(userPos(j,2) - rec_xyz_actual(2));
            err_z(j) = abs(userPos(j,3) - rec_xyz_actual(3));
            mean_err(j) = (err_x(j) + err_y(j) + err_x(j)) / 3;
        end
        
        subplot(2,2,1);
        hold on;
        if i ~= nargin-2
            plot(epochs-epochs(1), mean_err);
        else 
            plot(epochs-epochs(1), mean_err, "cyan");
        end
        xlabel("Time (s)");
        ylabel("Error (m)");
        xlim([0 epochs(nepoch)-epochs(1)]);
        grid("on");
        if i == nargin-2
            legend("LS", "Sigma-modeled WLS", "Exponential-modeled WLS", "Multilateration");
        end

        % title("Mean Position Error");
        % plot_settings(epochs, nepoch, nargin, i, method_name);
        
        subplot(2,2,2);
        hold on;
        plot(epochs-epochs(1), err_x);
        title("X axis Error");
        plot_settings(epochs, nepoch, nargin, i, method_name);
        
        subplot(2,2,3);
        hold on;
        plot(epochs-epochs(1), err_y);
        title("Y axis Error");
        plot_settings(epochs, nepoch, nargin, i, method_name);
        
        subplot(2,2,4);
        hold on;
        plot(epochs-epochs(1), err_z);
        title("Z axis Error");
        plot_settings(epochs, nepoch, nargin, i, method_name);

        fprintf('%s mean error is %d m\n',inputname(2+i), mean(mean_err));
        fprintf('%s X axis error is %d m\n',inputname(2+i), mean(err_x));
        fprintf('%s Y axis error is %d m\n',inputname(2+i), mean(err_y));
        fprintf('%s Z axis error is %d m\n',inputname(2+i), mean(err_z));
        fprintf('\n');
    end

    function plot_settings(epochs, nepoch, nargin, count_method, method_name)
        xlabel("Time (s)");
        ylabel("Error (m)");
        xlim([0 epochs(nepoch)-epochs(1)]);
        grid("on");
        
        if count_method == nargin-2
           switch (nargin-2)
                case 1
                    legend(method_name(1,:));
                case 2
                    legend(method_name(1,:), method_name(2,:));
                case 3
                    legend(method_name(1,:), method_name(2,:), ...
                        method_name(3,:));
                case 4
                    legend(method_name(1,:), method_name(2,:), ...
                        method_name(3,:), method_name(4,:));
                case 5
                    legend(method_name(1,:), method_name(2,:), ...
                        method_name(3,:), method_name(4,:), ...
                        method_name(5,:));
            end
        end
    end
end
