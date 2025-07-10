function retS = spec_ret_marginal(Mwin, noise, num_spec,npad) 
   
 
    Mw = quickscale(Mwin);
    N = length(Mw);
    Mw = Mw(:).';
    
    if noise        
        peri = 5;
        Mw = side_noise_1d(Mw, peri);
 
        Mw(real(Mw)<0) = 1e-6; 
 
        Mw = smooth(Mw,'sgolay');
        Mw = Mw(:).';
  
        Mw = super_gaussian_1d(Mw, .9,10);
 
        Mw(real(Mw)<0) = 1e-6; 
        
    end
    
    
    Mw = quickscale(Mw);
    
    coeff_rmv = .2;
  
    St2 = ifftc(1/2/pi*Mw,npad);
 
    St = sqrt(St2);
 
     
    imagSt = imag(St);          realSt = real(St);      % orginal real and imaginary parts of s(t)
   
    imagStc = imagSt;           realStc = realSt;       % corrected real and imaginary parts of s(t)
    
    pos_neg = ones(size(imagSt));                       % postive or negative coefficent for each compponent of St
      
    
    rmsl = [];   rmsr = [];   
 
    %%% %%% %%%
 
    c = [.09, .425, 1];
    
%%    
    
    [~, minima_ind] = findpeaks(-realSt);
    [~, maxima_ind] = findpeaks(realSt); 
    
  
    minima_ind(minima_ind < round(N * coeff_rmv)) = [];
    minima_ind(minima_ind > N - round(N * coeff_rmv)) = [];
    
    
    maxima_ind(maxima_ind < round(N * coeff_rmv)) = [];
    maxima_ind(maxima_ind > N - round(N * coeff_rmv)) = [];
    
    nix = length(maxima_ind);
    nim = length(minima_ind);
    
    minima_ind_org = minima_ind;
    maxima_ind_org = maxima_ind;
    
    %% check needed
    % making them in the maxima and minima in the same order 
    [~, indp1] =  max(abs(realSt));
    locus = find( minima_ind <  indp1 );
   
 
    % ordering the minima peakes based on the maxima intensity to their
    % right/ not a necessary step
    if length(locus) > 0 
        differ = locus(end) - find(maxima_ind == indp1);   
        minima_ind = circshift(minima_ind, -differ);   
    else
        differ = 0;
    end
    
    if differ < 0
        minima_ind(1:abs(differ)) = [];
        maxima_ind(1:abs(differ)) = [];
    elseif differ > 0
        minima_ind(end-abs(differ):end) = [];
        maxima_ind(end-abs(differ):end) = [];
    end
    % %
    
   if nix > nim
        maxima_ind(length(minima_ind)+1:end)=[];
   else
        minima_ind(length(maxima_ind)+1:end)=[];
   end    
    
    
    fR = abs(realSt);
    fR(indp1+1:end) = 0;
    min_dum =[];
    min_dum2 =[];
    min_dum3 = [];
    
    for i_peak = 1: min(log2(N) ,sum(maxima_ind<N/2+2))
 
        [~, indmax(i_peak)] = max(abs(fR(maxima_ind)));
    
        min_ind = minima_ind(minima_ind < maxima_ind(indmax(i_peak)));  
        min_dum = min_ind(end);
        
        min_ind2 = minima_ind(minima_ind > maxima_ind(indmax(i_peak)));  
        if ~isempty(min_ind2)
            min_dum2 = min_ind2(1);
        end
        
        min_dum3 = [min_dum3, min_dum];
        fR(maxima_ind(indmax(i_peak))) = 0;
        
    end
 
 %%
 
    minima_ind = min_dum3(min_dum3<N/2+1);
 
    % removing the points with large imaginary parts 
 
    
    if length(minima_ind)> 10
        minima_ind = minima_ind(1: round( 10));
    end   
    
  
    
    min_ind_left = minima_ind(minima_ind<N/2+1);
  
    
%% % right check
 
    minima_ind = minima_ind_org;
    maxima_ind = maxima_ind_org;
    
 
    % making them in the maxima and minima in the same order 
    [~, indp1] =  max(abs(realSt));
    locus = find( minima_ind <  indp1 );
   
    
    % ordering the minima peakes based on the maxima intensity to their
    % right/ not a necessary step
    if length(locus) > 0
        
        differ = locus(end) - find(maxima_ind == indp1);
    
        minima_ind = circshift(minima_ind, -differ);   
    else
        differ = 0;
    end
 
    
    if differ < 0
        minima_ind(1:abs(differ)) = [];
        maxima_ind(1:abs(differ)) = [];
    elseif differ > 0
        minima_ind(end-abs(differ):end) = [];
        maxima_ind(end-abs(differ):end) = [];
    end
    % %
 
   if nix > nim
        maxima_ind(length(minima_ind)+1:end)=[];
   else
        minima_ind(length(maxima_ind)+1:end)=[];
   end    
    
    
    fR = abs(realSt);
    fR(1:indp1) = 0;
    min_dum =[];
    min_dum2 =[];
    min_dum3 =[];
 
    for i_peak = 1: min(log2(N) ,sum(maxima_ind>N/2+1))
 
        [~, indmax(i_peak)] = max(abs(fR(maxima_ind)));
    
        min_ind = minima_ind( minima_ind > maxima_ind(indmax(i_peak)));  
                
        if ~isempty(min_ind)
            min_dum = min_ind(1) ;
        end
        
    
        min_ind2 = minima_ind( minima_ind < maxima_ind(indmax(i_peak)));  
        
        if ~isempty(min_ind2)
            min_dum2 = min_ind2(end) ;
        end
        
        
        fR(maxima_ind(indmax(i_peak))) = 0;
        
        min_dum3 = [ min_dum3,min_dum2];
   
    end
 
    %% right cont
    minima_ind = min_dum3(min_dum3>N/2+1);
    
    if length(minima_ind)> 11
        minima_ind = minima_ind(1: 1: 11 );
    end      
    
%     
    nn = 1;
    if length(minima_ind) > 10 
 
        fakeImag1 = abs(imagStc(minima_ind)); 
        fakeImag2 = abs(imagStc(minima_ind-1));
        fakeImag3 = abs(imagStc(minima_ind+1));
 
        if length(minima_ind)<3    
            remover = fakeImag1 > .85 * max(abs(imagStc))  & fakeImag2 > .98 * max(abs(imagStc)) & fakeImag3 > .98 * max(abs(imagStc));
        else
            remover = fakeImag1 > .35^nn * max(abs(imagStc))  & fakeImag2 > .35^nn * max(abs(imagStc));
        end
 
        [val, ind_dum] = sort(minima_ind);
        minima_ind( ind_dum(remover) ) = [];  
        nn = nn+.05;
    end
    
    
    if length(minima_ind)> 10%log2(N)-2
        minima_ind = minima_ind(1: 1: 10 );
    end      

       
    % removing the points with large imaginary parts 
 
     min_ind_right = minima_ind(minima_ind>N/2+1);
   
     %
 
%%
    
    add_l = length(min_ind_left);
    add_r = length(min_ind_right);
    
    num_left = 2^add_l;
    num_right = 2^add_r;
    
    
    pos_neg_left = ones(num_left,N);
    pos_neg_right = ones(num_right,N);
    
    
    ind_dum = sort(min_ind_left);
    for ii = 1:add_l
        ind0 = 1:2^add_l / 2^ii;
        ind1 = ind0;
        ind0=[];
        for jj = 1:2^(ii-1)
            ind0 = [ind0 , ind1];
            ind1 = ind1 + 2^add_l/2^(ii-1);
        end
 
        
        if ii < add_l
                pos_neg_left(ind0,ind_dum(ii):ind_dum(ii+1)-1) = ...
                    -pos_neg_left(ind0,ind_dum(ii):ind_dum(ii+1)-1);
        else
                pos_neg_left(ind0,ind_dum(ii):N/2+1) =  -pos_neg_left(ind0,ind_dum(ii):N/2+1);
        end
 
    end
 
    pos_neg_left(:,N/2+2:end) = fliplr(pos_neg_left(:,2:N/2));
    
  
    
    % % 
    
    ind_dum = sort(min_ind_right);
    for ii = 1:add_r
        ind0 = 1:2^add_r / 2^ii;
        ind1 = ind0;
        ind0=[];
        for jj = 1:2^(ii-1)
            ind0 = [ind0 , ind1];
            ind1 = ind1 + 2^add_r/2^(ii-1);
        end
 
       
        if ii < add_r
                pos_neg_right(ind0, ind_dum(ii):ind_dum(ii+1)-1) =  ...
                    -pos_neg_right(ind0, ind_dum(ii):ind_dum(ii+1)-1);
        else
                pos_neg_right(ind0,ind_dum(ii):N) =  -pos_neg_right(ind0,ind_dum(ii):N);
        end
 
    end
 
    pos_neg_right(:,2:N/2)  = fliplr(pos_neg_right(:,N/2+2:end));
    
   
    
     %% 
    
  
    for i_step =  5: N
 
         dum_realStc = realStc(1:i_step);
         dum_imagStc = imagStc(1:i_step);
 
         dum_realStc(i_step) = -realStc(i_step);
         dum_imagStc(i_step) = -imagStc(i_step);
         
         
         D3realStc = diff(diff(diff(realStc(1:i_step))));
         D3imagStc = diff(diff(diff(imagStc(1:i_step))));
 
 
         dum_D3realStc = diff(diff(diff(dum_realStc)));
         dum_D3imagStc = diff(diff(diff(dum_imagStc)));
         
         D12r = realStc(i_step) - realStc(i_step-1);
         D01r = realStc(i_step-1) - realStc(i_step-2);
         
         dum_D12r = dum_realStc(i_step) - realStc(i_step-1);
         dum_D01r = realStc(i_step-1) - realStc(i_step-2);
         
         D12im = imagStc(i_step) - imagStc(i_step-1);
         D01im = imagStc(i_step-1) - imagStc(i_step-2);
         
         dum_D12im = dum_imagStc(i_step) - imagStc(i_step-1);
         dum_D01im = imagStc(i_step-1) - imagStc(i_step-2);
         
 
         ci(1) = c(1);   cr(1) = c(1);
         ci(2) = c(2);   cr(2) = c(2);
         ci(3) = c(3);   cr(3) = c(3);
 
         
    
         wp(1) =  abs(D12r)^2 * cr(1) + ...
                  abs(D12r - D01r)^2 * cr(2)  + ...
                  abs(D3realStc(i_step-3))^2 * cr(3)  ;
         
         wn(1) =  abs(dum_D12r)^2 * cr(1) + ...
                  abs(dum_D12r - dum_D01r)^2 * cr(2) + ...
                  abs(dum_D3realStc(i_step-3))^2 * cr(3) ; 
                   
         wp(2) =  abs(D12im)^2 * ci(1) + ...
                  abs(D12im-D01im)^2 * ci(2) + ...
                  abs(D3imagStc(i_step-3))^2 * ci(3) ;
         
         wn(2) =  abs(dum_D12im).^2 * ci(1) + ...
                  abs(dum_D12im - dum_D01im)^2 * ci(2) + ...
                  abs(dum_D3imagStc(i_step-3))^2 * ci(3) ;
 
 
         if sum(wp) > sum(wn)
                realStc(i_step:end) = -realStc(i_step:end);
                imagStc(i_step:end) = -imagStc(i_step:end);
                pos_neg(i_step:end) = -pos_neg(i_step:end);        
         end
 
    end
    
    pos_neg1 = pos_neg;     pos_neg1(N/2+2:end) = fliplr(pos_neg1(2:N/2));
    pos_neg2 = pos_neg;     pos_neg2(2:N/2)= fliplr(pos_neg2(N/2+2:end));
 
    
    
    [Swc(1,:), rms0(1), rms0p(1)] = left_side(Mw, St, pos_neg1);
 
    [Swc(2,:), rms0(2), rms0p(2)] = right_side(Mw, St, pos_neg2);
    
 
% % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA         
% 
% 
        for i_peak = num_left:-1:1
            pos_neg11 = pos_neg1.*pos_neg_left(i_peak,:);
            [Swcl(i_peak, :),rmsl(i_peak), rmslp(i_peak)] = left_side(Mw, St, pos_neg11);            
        end
        
        for i_peak = num_right:-1:1
            pos_neg22 = pos_neg2.*pos_neg_right(i_peak,:);
            [Swcr(i_peak, :),rmsr(i_peak), rmsrp(i_peak)] = right_side(Mw, St, pos_neg22); 
 
        end
        
                
        rms_total = [ rmsr, rmsl];
        rms_total_p = [rmsrp, rmslp];
        S_total = [ Swcr; Swcl];
   
 
    %%             
                              
            [condf, order] = sort(rms_total);
            num_s = length(order);
 
            
            V = 0;
            if num_s > 2     
                V = 3;
            end    
 
            num_asg = 1;
            
                if  num_s <= 4
                   
                    Swf = S_total;
 
                elseif  num_s <= 12
                        ind = min( round(num_s)); 
                        Swf(1,:) = quickscale(S_total(order(1),:));
                        S_total_p = S_total(order(2:ind),:);
                        order = order(2:ind);
                        
 
 
                        for ii = 1:V               
                            num_asg = num_asg+1; 
                            delta = ones(1,length(order));
                            for jj = 1: length(order)
 
                                    for kk = 1:num_asg-1
                                         delta_2 = mean(abs(quickscale(S_total_p(jj,:)) - quickscale(Swf(kk,:))).^2); 
                                         delta(jj) = delta(jj) *delta_2;
                                    end
 
                               
                            end
                            [~,ind] = sort(delta);
                            Swf(ii+1,:) = quickscale(S_total_p(ind(end),:));        %selecting the most diffirent ones
                            order(ind(end))=[];
                            S_total_p(ind(end),:)=[];
                        end
               
                else
                    if num_s < 16
                        coeff_up = 12/num_s;
                    elseif  num_s < 26
                        coeff_up = 12/num_s;     
                    elseif num_s < 40
                        coeff_up = 15/num_s;   
                    elseif num_s < 63
                        coeff_up = 18/num_s;                                      
                    else
                        coeff_up = 25/num_s;
                    end
                        
                    [~, order_p] = sort(rms_total_p);
                    
                    
                    
                    diff12 = sqrt(mean(abs(S_total(order(1),:) - S_total(order(2),:)).^2));

                    if diff12>.08
                        
                        
                        
                        [order_union, io, iop] = intersect(order(3:round(coeff_up*num_s)), order_p(1:round(coeff_up*num_s)));
                        
                        while length(io) < 8
                            coeff_up = coeff_up*1.1;
                            [order_union, io, iop] = intersect(order(3:round(coeff_up*num_s)), order_p(1:round(coeff_up*num_s)));
                        end
                        
                        Swf(1:2,:) = S_total(order(1:2),:);
                        S_total_u = S_total(order(2+io),:);
                        [~, order_n] = sort(rms_total( order(2+io)));
                        num_asg = 2;
                        start = 2;

                        
                    else    
                        
                        [order_union, io, iop] = intersect(order(2:round(coeff_up*num_s)), order_p(1:round(coeff_up*num_s)));
                        
                        while length(io) < 8
                            coeff_up = coeff_up*1.1;
                            [order_union, io, iop] = intersect(order(2:round(coeff_up*num_s)), order_p(1:round(coeff_up*num_s)));
                        end
                            
                        
                        Swf(1,:) = S_total(order(1),:);
                        S_total_u = S_total(order(1+io),:);
                        [~, order_n] = sort(rms_total(order(1+io)));
                        num_asg = 1;
                        start = 1;
                    end
                  
                 
                    S_total_p = S_total_u(order_n,:);
 
 
                    for ii = start:V                 
                        delta = ones(1,length(order_n));
 
                        for jj = 1: length(order_n)
 
                                for kk = 1:num_asg
                                     delta_2 = mean(abs(quickscale(S_total_p(jj,:)) - quickscale(Swf(kk,:))).^2); 
                                     delta(jj) = delta(jj)*sqrt(delta_2);
                                end
 
                            
                        end
                        [~,ind] = sort(delta);
                         
                        
                        Swf(ii+1,:) = quickscale(S_total_p(ind(end),:));        %selecting the most diffirent ones
                        order_n(ind(end))=[];
                        S_total_p(ind(end),:) = [];
                        num_asg = num_asg+1; 
 
                    end
         
           
                end
 
    
 
    
    if num_spec == 1
        retS = abs(Swf(1,:));
    else 
        retS = abs(Swf);
    end
       
    if   size(retS,1) == 3
        retS = [abs(Swf); abs(Swf(1,:))];
        
    end
    
 
end
 
 

