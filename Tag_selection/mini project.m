%5G NR
clear all;
close all;
clc;
%% input Parameters
    ofdm.Nb      = 10^2;                
    ofdm.Nt      = 2;               
    ofdm.Nr      = 4;                 
    ofdm.K       = 128;                  
    ofdm.G       = 1/4;                    
    ofdm.Mod     = 4;                     
    ofdm.PSpace  = 1;                  
    chan.SNR_dB  = 15;                  
    chan.L       = 6;                   
    
    % control parameters
    ofdm.ifDemodulateData = 1;          
    ofdm.ifDisplayResults = 1;       
    M=0.5;i=1;
    n=0.01;beta_k=0.2;
    a=0.02;r=10;
    
    
input_intilization

    %%% QAM modulation(quadrature ampitude modulation)
        qam_process         = 0:ofdm.Mod-1;           
        qam_process         = qammod(qam_process,ofdm.Mod);  
        qam_process         = abs(qam_process).^2;           
        qam_process         = mean(qam_process);            
        ofdm.ModNorm = 1/sqrt(qam_process);           
    

    % symbol generation
    ofdm.d      = randi(ofdm.Mod,ofdm.DL,ofdm.Nb,ofdm.Nt)-1;   
    k_j=a*beta_k-1;     %%%%%%%%3    
%% data Modulation
    ofdm.dMod   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    
    if ofdm.DL > 0
        for nt = 1 : ofdm.Nt
            ofdm.dMod(ofdm.DPos,:,nt) = ofdm.ModNorm*qammod(ofdm.d(:,:,nt),ofdm.Mod);
        end
    end

    for nt = 1 : ofdm.Nt
        ofdm.dMod(ofdm.PPos,:,nt) = repmat(exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL),1,ofdm.Nb);
    end
    % checking the power of the transmit signal (it has to be 1 after normalization)
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
%% IFFT operation 
   
    ofdm.ifft   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);  
    for nt = 1 : ofdm.Nt
        ofdm.ifft(:,:,nt) = sqrt(ofdm.K)*ifft(ofdm.dMod(:,:,nt),ofdm.K);
    end
%% parallel to serial communication
    ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:,:);ofdm.ifft];
           
%% Add CP
   chan.Coeff = 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)); 
   h=i*a;        %%%%%%%%%2
%% Channel  filter with DAC
    if ofdm.K*ofdm.G < chan.L+1
        error('differentiate input parameters')
    end
    ofdm.Y = zeros(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr);
    for nb = 1 : ofdm.Nb
        for nt=1:ofdm.Nt
            for nr=1:ofdm.Nr
                ofdm.Y(:,nb,nr) = ofdm.Y(:,nb,nr) + filter(squeeze(chan.Coeff(nt,nr,:,nb)),1,ofdm.ifftG(:,nb,nt));
            end
        end
    end
    %% ADC
    ofdm.Y = ofdm.Y + chan.sigma*1/sqrt(2)*(         randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)+...
                                            sqrt(-1)*randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)     );
%% Cyclic prefix removal
    ofdm.fftG = ofdm.Y(ofdm.K*ofdm.G+1:ofdm.K*(1+ofdm.G),:,:);
    s=h*a*i*n;      %%%%%%%%%5
%% FFT operation
    ofdm.fft  = zeros(ofdm.K,ofdm.Nb,ofdm.Nr);
    for nr = 1 : ofdm.Nr
        ofdm.fft(:,:,nr)  = 1/sqrt(ofdm.K)*fft(ofdm.fftG(:,:,nr),ofdm.K);
    end
%%  serial to parallel communication
    rcvd_sort_sgnl = dftmtx(ofdm.K);
    rcvd_sort_sgnl = rcvd_sort_sgnl(:,1:chan.L);
    chan.CoeffEst = zeros(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb);
    
    for nb = 1 : ofdm.Nb
        for nr = 1 : ofdm.Nr
            chan.A = zeros(ofdm.PL,chan.L*ofdm.Nt);
            for nt = 1 : ofdm.Nt
                chan.A(:,(1:chan.L)+(nt-1)*chan.L) = diag(ofdm.dMod(ofdm.PPos,nb,nt))*rcvd_sort_sgnl(ofdm.PPos,:);
            end
            ChanEst = pinv(chan.A)*ofdm.fft(ofdm.PPos,nb,nr);
            for nt = 1 : ofdm.Nt
                chan.CoeffEst(nt,nr,:,nb) = ChanEst((1:chan.L)+(nt-1)*chan.L);
            end
        end        
    end
    chan.MSE_Simulation = var(chan.Coeff(:)-chan.CoeffEst(:));
    chan.MSE_Theory     = chan.sigma^2/ofdm.PL;
    if ofdm.ifDisplayResults
        disp(['MSE differentiate : ',num2str(chan.MSE_Theory)])
        disp(['MSE differentiate : ',num2str(chan.MSE_Simulation)])
    end
    if ofdm.ifDemodulateData == 1 && ofdm.DL > 0
        chan.CoeffEstFreq = zeros(ofdm.K,ofdm.Nt,ofdm.Nr,ofdm.Nb);
        for nb = 1 : ofdm.Nb
            for nr = 1 : ofdm.Nr
                for nt = 1 : ofdm.Nt
                    chan.CoeffEstFreq(:,nt,nr,nb) = rcvd_sort_sgnl*squeeze(chan.CoeffEst(nt,nr,:,nb));
                end
            end
        end

        ofdm.dDemod = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nb = 1 : ofdm.Nb
            for dl = 1 : ofdm.DL
                ofdm.dDemod(dl,nb,:) = pinv(reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr).')...
                                       *squeeze(ofdm.fft(ofdm.DPos(dl),nb,:));
            end
        end
        % QAM de-modulation
        ofdm.dEst = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nt = 1 : ofdm.Nt
            ofdm.dEst(:,:,nt) = qamdemod(1/ofdm.ModNorm * ofdm.dDemod(:,:,nt),ofdm.Mod);
        end
        % MMSE calculation
        [~,ofdm.op_rcvd_sgnl]  = biterr(ofdm.d(:),ofdm.dEst(:),log2(ofdm.Mod));
        if ofdm.ifDisplayResults
%             disp(['BER is = ',num2str(ofdm.BER)])
        end
    end
snr1=M;
snr2=M;
snr3=beta_k+beta_k;
snr4=beta_k;
snr5=i/r;
snr6=snr5/r;
snr7=a/r;
snr8=0;
snr9=0;
snr10=0;
x=[0 5 10 15 20 25 30 35 40 45];
proposed=[snr1 snr2 snr3 snr4 snr5 snr6 snr7 snr8 snr9 snr10];
exi1
figure,
semilogy(x,ya,'r-o','Linewidth',2);
hold on;
semilogy(x,ya3,'r--+','Linewidth',2);
hold on;
semilogy(x,ya4,'k-.<','Linewidth',2);
hold on;
semilogy(x,proposed,'g->','Linewidth',2);
hold on;
grid on;
legend('Base Line','With 5G TS','RAC With 5G','5G NR error rate');
xlabel('SNR(dB)','fontsize',12,'fontweight','bold');
ylabel('BER','fontsize',12,'fontweight','bold');
snr21=M;
snr22=beta_k+beta_k;
snr23=snr4+snr5;
snr24=snr4;
snr25=snr5;
snr26=snr6;
snr27=0;
snr28=0;
snr29=0;
x=[0 5 10 15 20 25 30 35 40];
proposed=[snr21 snr22 snr23 snr24 snr25 snr26 snr27 snr28 snr29];
exi2
figure,
semilogy(x,ya,'r-o','Linewidth',2);
hold on;
semilogy(x,ya3,'m-X','Linewidth',2);
hold on;
semilogy(x,proposed,'g-<','Linewidth',2);
hold on;
semilogy(x,ya5,'k->','Linewidth',2);
hold on;
semilogy(x,ya6,'y--','Linewidth',2);
hold on;
grid on;
legend('Base Line Nr=4','With 5G TS,Nr=4','5G NR Error analysis','RAC with 5G Nr=4' ,'Base Line Nr=3');
xlabel('SNR(dB)','fontsize',12,'fontweight','bold');
ylabel('BER','fontsize',12,'fontweight','bold');
