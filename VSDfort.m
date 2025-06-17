function [Utmp,d,err] = VSDfort(T,R,S)
% Coded by Panqi Chen, 2023.
%   CPD4v Vandermonde-constrained CPD using spatial smoothing，three of the
%   factor matrices(here we specify 1-th,2-th and 3-th are vandermonde) of the CPD are Vandermonde-constrained;
%   [U,d,fit_error] = cpd3v_ss2(T,R,S) computes the factor matrices U{1}, U{2}, 
%   and U{3} belonging to a canonical polaydic decomposition of a third-
%   order tensor T with rank R of which the 1,2-th factor matrix is 
%   Vandermonde-constrained. The decomposition is computed via an 
%   eigenvalue decomposition.
%    
%   ~~~~~~~~~~~~~~~~~~~~output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   U returns the estimated factor matrix.
%   d returns the estimated Vandermonde generators.
%   error returns the tensor reconstructed error.

    % Check inputs. T=X  
    if getorder(T) ~= 4
        error('cpd4v:T','Only fourth order tensors are allowed.');    
    end
    size_tens=size(T);% S=[12,36,2]
    K1=S(1);K2=S(2);K3=S(3);
    L1=size_tens(1)+1-K1;
    L2=size_tens(2)+1-K2;
    L3=size_tens(3)+1-K3;
  

  % Apply spatial smoothing (in structured mode).
    Y = hankelize(T,'Dim',1,'Ind', K1); 
    YY = hankelize(Y,'Dim',3,'Ind', K2); 
    YYY= hankelize(YY,'Dim',5,'Ind', K3); 
     Ymat2=tens2mat(YYY,[5 3 1],[7 6 4 2]);   %norm(Ymat2-X_mat,2)
      
    % Compute the first group of vandermonde factors
      
      Utmp = cell(1,4);
      [Uy,~,~] = svd(Ymat2);  %svd
      Zhat = Uy(1:end-K2*K3,1:R)\Uy(K2*K3+1:end,1:R); 
      [M1,D1] = eig(Zhat);
      d1 = diag(D1);
      d1=d1./abs(d1); 
      Utmp{1} = struct_vander(d1,[],size_tens(1)-1,true);
        
        
        
        tmp = bsxfun(@times, reshape(Uy(:,1:R)*M1,[K2*K3 K1 R]), reshape(conj(Utmp{1}(1:K1,:)),[1 K1 R]));  
        Utmp2_K2 = squeeze(sum(tmp,2));
        
        Zhat2 =[];
        for i=1:R
            ztmp=Utmp2_K2(1:end-K3,i)\ Utmp2_K2(K3+1:end,i); 
            Zhat2 =[Zhat2;ztmp];
        end  
         d2=Zhat2./abs(Zhat2); %% 归一化
         Utmp{2} = struct_vander(d2,[],size_tens(2)-1,true);
         
         
         
        tmp2 = bsxfun(@times, reshape(Utmp2_K2(:,1:R),[K3 K2 R]), reshape(conj(Utmp{2}(1:K2,:)),[1 K2 R]));    
        Utmp3_K3 = squeeze(sum(tmp2,2));
         %determine vande generator 2
        d3=[];
        for j=1:R
             dtemp=Utmp3_K3(1:end-1,j)\Utmp3_K3(2:end,j);
             d3=[d3;dtemp];
        end       
       d3=d3./abs(d3); 
        Utmp{3} = struct_vander(d3,[],size_tens(3)-1,true);
   
        
        UHU = (Utmp{1}'*Utmp{1}).*(Utmp{2}'*Utmp{2}).*(Utmp{3}'*Utmp{3});
        UHUinv = pinv(UHU);
        Utmp{4} = mtkrprod(T,Utmp,4)*conj(UHUinv);      
        d=[d1,d2,d3];
        Y_recon=cpdgen(Utmp);
        err=norm(Y_recon(:)-T(:),2);
        
end