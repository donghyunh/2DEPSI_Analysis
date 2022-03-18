
function  [MetImageC_TD, MetImageC_TD_ZF, MetImageC_FD_ZF] = Zero_filling(MetImageC,zf_factor,starting_pt,ending_pt);

[a,b,c,d] = size(MetImageC);
MetImageC_TD= zeros(a,b,c,d);
MetImageC_FD= zeros(a,b,c,d);

for x = 1:a;
    for y = 1:b;
        for z = 1:d;
        MetImageC_TD(x,y,:,z) = ifftshift(ifft(ifftshift(MetImageC(x,y,:,z)))); %1*1*256
        MetImageC_TD_ZF(x,y,:,z) = cat(3,MetImageC_TD(x,y,:,z), zeros(1,1,zf_factor));
        MetImageC_FD_ZF(x,y,:,z) = fftshift(fft(fftshift(MetImageC_TD_ZF(x,y,:,z))));
        
        end
    end
end

% size(MetImageC_TD)
% size(MetImageC_TD_ZF)
% size(MetImageC_FD_ZF)
% 
% figure,
% subplot(2,1,1)
% plot(squeeze(squeeze(real(MetImageC_TD_ZF(5,3,:,4))))); hold on
% plot(squeeze(squeeze(imag(MetImageC_TD_ZF(5,3,:,4))))); 
% subplot(2,1,2)
% plot(squeeze(squeeze(abs(MetImageC_TD_ZF(5,3,:,4)))));
% 
% figure,
% subplot(2,1,1)
% plot(squeeze(squeeze(real(MetImageC_FD_ZF(5,3,:,4))))); hold on
% plot(squeeze(squeeze(imag(MetImageC_FD_ZF(5,3,:,4)))));
% subplot(2,1,2)
% plot(squeeze(squeeze(abs(MetImageC_FD_ZF(5,3,:,4)))));

end
