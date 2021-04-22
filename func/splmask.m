%%=========     SCALING FACTOR CONSIDERATION    ============%%
%its purpose is to energy preserving by L2-norm coefficients
%that is: normalization should be the one satisfying Bessel inequality

%====   the downsampled version ============================%% 
%the normalization factor is 1/sqrt(2)
%with respect to the mask in the refinement or wavelet equation

%====   the undecimated version ============================%% 
%the normalization factor is 1/2
%with respect to the mask in the refinement or wavelet equation

%%%=========================================================%%
%please refer to p149 <<A wavelet tour of signal processing>>

%the mask for linear or cubic spline wavelets, and Harr wavelet
function m = splmask(order, level, version)

pad=zeros(1,2^(level-1)-1);
if order == 0
    %linear spline frame masks
    if version == 'U' %undecimated
        phi=1/2*[zeros(1,level==1)  1 pad 1];       
        psi=1/2*[zeros(1,level==1) 1 pad -1]; 
    else  %version == 'D'      
        phi=1/2*[1 pad 1 zeros(1,level==1)];
        psi=1/2*[-1 pad 1 zeros(1,level==1)]; 
    end
    m=[phi;psi];
elseif order==1    
    %linear spline frame masks
    if version == 'U' %undecimated
        phi=1/4*[1 pad 2 pad 1];       %refinement mask
        psi1=sqrt(2)/4*[1 pad 0 pad -1];    %first wavelet mask
        psi2=1/4*[-1 pad 2 pad -1];         %second wavelet mask
    else  %version == 'D'      
        phi=sqrt(2)/4*[1 pad 2 pad 1];     
        psi1=1/2*[1 pad 0 pad -1];         
        psi2=sqrt(2)/4*[-1 pad 2 pad -1];  
    end
    m=[phi;psi1;psi2];    
elseif order==3    
    %cubic spline frame masks
    if version == 'U' %undecimated
        phi=[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];            %refinement mask
        psi1=[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];             %first wavelet mask
        psi2=sqrt(3)/sqrt(2)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];%second wavelet mask
        psi3=[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];             %third wavelet mask
        psi4=[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];         %fourth wavelet mask
    else
        phi=sqrt(2)*[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];            
        psi1=sqrt(2)*[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];             
        psi2=sqrt(3)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];                
        psi3=sqrt(2)*[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];             
        psi4=sqrt(2)*[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];         
    end
    m=[phi;psi1;psi2;psi3;psi4];
else
    error('wrong spline order');
end
end

%% Original Version which has bugs on Haar
% function m = splmask(order, level, version)
% 
% pad=zeros(1,2^(level-1)-1);
% if order == 0
%     %linear spline frame masks
%     if level == 1
%         if version == 'U' %undecimated
%             phi=1/2*[0 1 pad 1];       
%             psi=1/2*[0 1 pad -1]; 
%         else  %version == 'D'      
%             phi=1/2*[1 pad 1 0];
%             psi=1/2*[-1 pad 1 0]; 
%         end
%     else
%         if version == 'U' %undecimated
%             phi=1/2*[1 pad 1];       
%             psi=1/2*[1 pad -1]; 
%         else  %version == 'D'      
%             phi=1/2*[1 pad 1];
%             psi=1/2*[-1 pad 1]; 
%         end        
%     end
%     m=[phi;psi];
% elseif order==1    
%     %linear spline frame masks
%     if version == 'U' %undecimated
%         phi=1/4*[1 pad 2 pad 1];       %refinement mask
%         psi1=sqrt(2)/4*[1 pad 0 pad -1];    %first wavelet mask
%         psi2=1/4*[-1 pad 2 pad -1];         %second wavelet mask
%     else  %version == 'D'      
%         phi=sqrt(2)/4*[1 pad 2 pad 1];     
%         psi1=1/2*[1 pad 0 pad -1];         
%         psi2=sqrt(2)/4*[-1 pad 2 pad -1];  
%     end
%     m=[phi;psi1;psi2];    
% elseif order==3    
%     %cubic spline frame masks
%     if version == 'U' %undecimated
%         phi=[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];            %refinement mask
%         psi1=[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];             %first wavelet mask
%         psi2=sqrt(3)/sqrt(2)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];%second wavelet mask
%         psi3=[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];             %third wavelet mask
%         psi4=[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];         %fourth wavelet mask
%     else
%         phi=sqrt(2)*[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];            
%         psi1=sqrt(2)*[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];             
%         psi2=sqrt(3)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];                
%         psi3=sqrt(2)*[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];             
%         psi4=sqrt(2)*[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];         
%     end
%     m=[phi;psi1;psi2;psi3;psi4];
% else
%     error('wrong spline order');
% end   
% end