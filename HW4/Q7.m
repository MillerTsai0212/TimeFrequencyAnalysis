clear all

t = [0: 0.01: 10];
x = 0.2*t + cos(2*pi*t) + 0.4*cos(10*pi*t);
thr = 0.2;
tic													% Timer Start
y = hht(x, t, thr);
toc													% Timer End



function y = hht(x, t, thr)
n = 1;
dt = diff(t(1:2));
max_iterations = 100;  % Set a maximum number of iterations to avoid infinite loop
iteration_count = 0;   % Initialize iteration counter
while (iteration_count < max_iterations)
    if length(findpeaks(x)) <= 3
        y(n,:) = x;
        break
    end
    
    temp = x;
    test = 1;
    k = 1;
    max_iterations_k = 3;  % Set a maximum number of iterations for k
    iteration_count_k = 0;   % Initialize iteration counter

    while (iteration_count_k < max_iterations_k)
        iteration_count_k = iteration_count_k + 1;
    
        [max, maxloc] = findpeaks(temp);         	% Step02
        peaks = interp1((maxloc-1)*dt, max, t, 'linear', 'extrap');  
        
        
        [neg_min, minloc] = findpeaks(temp*(-1));	% Step04
        min = -1 * neg_min;
        
        dips = interp1((minloc-1)*dt, min, t, 'linear', 'extrap');     	% Step05
        z = (peaks + dips) / 2;                    	% Step06_1
        h = (temp - z);								% Step06_2
        
        % Step07
        hpeaks = findpeaks(h);
        hdips = -findpeaks(-h);
        
        test = 0;
        for i = 1:(length(hpeaks)-1)
            if ((hpeaks(i) <= 0) || (hdips(i) >= 0) || (abs((hpeaks(i) + hdips(i))/2) >= thr))
                temp = h;
                test = 1;
                break
            end
        end
    end
    y (n,:) = h;
    % Step08
    x = x - h;
    n = n + 1;
end
plot(t, y)
end

