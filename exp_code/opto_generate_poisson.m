function stim_signal = opto_generate_poisson(f_sample, pulse_frequency, duration, amplitude, pulse_width, n_channels, refractory_period )
%   Detailed explanation goes here

if pulse_width>1-10^-10
    error('Wrong pulse width too big!')
end

%pulse width represents the bin size to determine which when/if pulse
%occurs
dt = pulse_width+refractory_period;
%calculate the number of bins to determine when pulses occur
sz = round(duration/dt);
%random numbers determine which bins will receive stimulation
x = rand(1,sz);
%random numbers below pulse_frequency*dt will have stim_pulse
pulse_t = (x <= pulse_frequency*dt);
time = [1:size(pulse_t,2)]*dt;

stim_signal = zeros(1,floor(size(pulse_t,2)*dt*f_sample));
stim_pulse = [amplitude*ones(1,ceil(pulse_width*f_sample)),zeros(1,ceil(refractory_period*f_sample))];
for i = 1:size(pulse_t,2)
    if pulse_t(i) == 0
        stim_signal(floor((i-1)*f_sample*dt+1):floor(i*f_sample*dt)) = 0;
    elseif pulse_t(i) == 1
        i1 = floor((i-1)*f_sample*dt+1);
        i2 = floor(i*f_sample*dt);
        stim_signal(i1:i2) = stim_pulse(1:i2-i1+1);
    end
end

% t   = 0:1/f_sample:duration;
% stim_signal = amplitude*interp1(time,double(pulse_t),t,'previous');
% dP = pulse_width/2:1/pulse_frequency:duration+pulse_width/2;
% yP = amplitude*pulstran(t,dP,'rectpuls',pulse_width);

% stim_signal = yP;
stim_signal = repmat(stim_signal,n_channels,1);

end
