clc;
clear;

// TASK 1

function [output] = dft_my(data)
    N = length(data);
    for k = 0:N-1
        s = 0.;
        for n = 0:N-1
            s = s + data(n+1) * exp(-2 * %i * %pi * n * k / N);
        end
        output(k+1) = s;
    end
endfunction

function [output] = fft_my(data)
    output = conj(fftReq(data)');
endfunction

function [output] = fft_req(data)
    N = length(data);
    if N == 1 then
        output = data;
    else
        data_even = fft_req(data(1:2:N));
        data_odd = fft_req(data(2:2:N));
        fact = exp(-2 * %i * %pi / N * (0:N-1));
        left = [data_even + fact(1:N/2) .* data_odd];
        right = [data_even + fact(N/2+1:N) .* data_odd];
        output = [left, right];
    end
endfunction


// TASKS 2-3

//fm = 10;
//fs = 125;
//m = 2.5; // without leakage = 2, with leakage = 2.5
//
//t = 0.0001 : 1/fs : m/fm;
//t_padded = 0.0001 : 1/fs : 1/fs * 64;
//
//x = cos(2 * %pi * fm/2 * t); // fm/2 for with leakage2
//
//N = (m * fs/fm);
//
//X1 = fft(x)(1:floor(N));
//spectral_leakage = abs(X1);
//zero_padded_x = resize_matrix(x, 1, 64) // zero padding
//spectral_interpolation = fft(zero_padded_x, 1) // ifft
//
//k = 0 : N-1;
//f = k * fs/N;
//
//f_norm = f ./ f(length(f))
//
//f = k / N;
//
//figure(1),
//subplot(221),
//plot2d(t, x),
//xlabel('time'),
//ylabel( 'x(n) '),
//title( 'With spectral leakage : m = 2, f=10 and Fs=125 Hz '),
//subplot(223),
//plot2d(f, spectral_leakage),
//xlabel('freq in Hz'), ylabel('Mag');
//subplot(222),
//plot2d(t_padded, abs(zero_padded_x)),
//subplot(224),
//plot(abs(spectral_interpolation))


fm1 = 190;
fs1 = 380;
m1 = 2;

t1 = 0 : 1/fs1 : 0.5;
x1 = 0.5 * cos(2 * %pi * fm1 * t1);

freq1 = t1 ./ t1(length(t1)) * fs1;
fft1 = fft(x1);

fm2 = 10;
fs2 = 200;
m2 = 2;

t2 = 0 : 1/fs2 : 0.5;
x2 = 2 * cos(2 * %pi * fm2 * t2);

freq2 = t2 ./ t2(length(t2)) * fs2;
fft2 = fft(x2);


figure(1),
subplot(221),
plot2d(t1, x1),
xlabel('time t'),
ylabel('cos(t)'),
title('T = 0.5, F = 190Hz, A = 0.5, Fs = 380 S/s'),
subplot(223),
plot2d(t2, x2),
xlabel('time t'),
ylabel('cos(t)');
title('T = 0.5, F = 10Hz, A = 2, Fs = 200 S/s');
subplot(222),
plot2d(freq1, abs(fft1)),
xlabel('frequency f'),
ylabel('magnitude(f)'),
title('fft from T = 0.5, F = 190Hz, A = 0.5, Fs = 380 S/s'),
subplot(224),
plot2d(freq2, abs(fft2)),
xlabel('frequency f'),
ylabel('magnitude(f)'),
title('fft from T = 0.5, F = 10Hz, A = 2, Fs = 200 S/s'),
