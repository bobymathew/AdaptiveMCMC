load('D:\patric\workshop_data.txt');
pedigree=workshop_data;
%y=y_a_d_e;
line_level =pedigree(:,1);
sire=pedigree(:,2);
dam=pedigree(:,3);
inbreed=pedigree(:,4);
y=pedigree(:,5);
%number of lines
n=size(pedigree,1);
records =n; 
A=zeros(n);
test=5000;
niter=50000;
burnin=2000;
div=50;
kA=1;lamA=1;
kD=1;lamD=1;
kE=1;lamE=1;
%programming the algorithm for A
for k=1:n
    for l=k:n
        if (k~=l)
            %both parents are unknown
            if sire(l)==0
                if dam(l)==0
                z=0;
                A(k,l)=z+A(k,l);
                A(l,k)=A(k,l);
                end
            end
            %one parent is known
            if sire(l)==0
                if dam(l)~=0
                z=0.5*A(k,dam(l));
                A(k,l)=z+A(k,l);
                A(l,k)=A(k,l);
                end
            end
            if sire(l)~=0
                if dam(l)==0
                z=0.5*A(k,sire(l));
                A(k,l)=z+A(k,l);
                A(l,k)=A(k,l);
                end
            end
            %both parents are known
            if sire(l)~=0
                if dam(l)~=0
                z=0.5*[A(k,sire(l))+A(k,dam(l))];
                A(k,l)=z+A(k,l);
                A(l,k)=A(k,l);
                end
            end
       
        else
            
              A(k,l)=1+inbreed(k);
            end
        end
    end



D=zeros(n);

%programming the algorithm for D
for b=1:n
    for c=1:n
        if (b~=c)
            %parents of both plants are unknown
            if sire(b)==0
            if dam(b)==0
            if sire(c)==0
            if dam(c)==0
                D(b,c)=0;
            end
            end
            end
            end
            
            %one parent of one plant is unknown
            if sire(b)==0
            if dam(b)~=0
            if sire(c)~=0
            if dam(c)~=0
                D(b,c)=0;
            end
            end
            end
            end
            
            if sire(b)~=0
            if dam(b)==0
            if sire(c)~=0
            if dam(c)~=0
                D(b,c)=0;
            end
            end
            end
            end
            
            if sire(b)~=0
            if dam(b)~=0
            if sire(c)==0
            if dam(c)~=0
                D(b,c)=0;
            end
            end
            end
            end
            
            if sire(b)~=0
            if dam(b)~=0
            if sire(c)~=0
            if dam(c)==0
                D(b,c)=0;
            end
            end
            end
            end
            
            %both parents of both plants are known
            if sire(b)~=0
            if dam(b)~=0
            if sire(c)~=0
            if dam(c)~=0
                D(b,c)=[[A(sire(b),sire(c))*A(dam(b),dam(c))]+[A(sire(b),dam(c))*A(dam(b),sire(c))]]/4;
            end
            end
            end
            end
        else
            D(b,c)=1;
        end
    end
end

X=(ones(records,1));  


vv=var(y);


y=(y-mean(y))/sqrt(var(y));

A_inv = inv(A);
A_sqrtm = sqrtm(A);
D_inv = inv(D);
D_sqrtm = sqrtm(D);
loc=0;
n_fix=1+loc;
Z=sparse(zeros(records,n));

for i=1:records
    for j=1:n
        if j==line_level(i)
            Z(i,j)=1;
        else
            Z(i,j)=0;
        end
    end
end


W=[X Z Z];
Wy=W'*y;
t=zeros(n_fix+n+n,1);
k_ratA=1;
k_ratE=1;

bl11=zeros(n_fix,n_fix);
bl12=zeros(n_fix,n);
bl21=zeros(n,n_fix);
bl13=zeros(n_fix,n);
bl31=zeros(n,n_fix);
bl23=zeros(n,n);
bl32=zeros(n,n);

E=[bl11 bl12 bl13; bl21 A_inv bl23; bl31 bl32 D_inv];

Cw=sparse(W'*W);
C=Cw+E;
fva=[];fve=[];fvp=[];ft=[];
fva_bn=[];fve_bn=[];fvp_bn=[];
fva_c=[];fve_c=[];fvp_c=[];
bn=1;
va=1;vp=1;ve=1;
si_a=1;si_p=1;si_e=1;
mn=0;
aa=0;
mmn=0;
acept=[];
count=0;
res=0;
flg=0;
kA=1;kD=1;kE=1;
kA=kA+n/2;kD=kD+n/2;kE=kE+records/2;%new K values for gamma
%starting the updating process (blocked gibbs sampler for Gaussian linear models)
for iteration=0:niter
    if(iteration<=test)
    if(mod(iteration,div)==0)
        a_star1= A_sqrtm*randn(n,1);         %generate a (first step)
        a_star=0+((1/sqrt((si_a)))*a_star1');%generate a (second step)  
        p_star1= D_sqrtm*randn(n,1);         %generate p (first step)
        p_star=0+((1/sqrt((si_p)))*p_star1'); 
        z_star=(Z*a_star')+(Z*p_star')+((1/sqrt((si_e)))*randn(records,1));
        Wf=W'*(y-z_star);
        LH=[zeros(1,n_fix), a_star, p_star]';
                
                    t=LH+(C\Wf);%inv(C)*Wf
                    
    else
                for i=1:n+n_fix+n
                    nsum=0;
                    for j=1:n+n_fix+n
                        nsum=nsum+(C(i,j)*t(j));
                    end
                    nsum=nsum-(C(i,i)*t(i));
                    nmean=(Wy(i)-nsum)/C(i,i);
                    vari=1/(si_e*C(i,i));
                    theta=nmean +(sqrt(vari)*randn(1,1));
                    t(i)=theta;
                end
                
    end
    
    
      
    lamA=0.001;lamD=0.001;lamE=0.001; %starting values as the variance of the phenotypic observations
    sa=t((n_fix+1):(n+n_fix))'*A_inv*t((n_fix+1):(n+n_fix));
    sp=t((n_fix+n+1):(n_fix+n+n))'*D_inv*t((n_fix+n+1):(n_fix+n+n));
    se=(y-W*t)'*(y-W*t);

    lamA=lamA+sa/2;lamD=lamD+sp/2;lamE=lamE+se/2;
    end

    if(iteration<=test) 
        
         ch_a=chi2rnd(2*kA);        %sampling additive variance
         ch_p=chi2rnd(2*kD);   %sampling interaction variance
         ch_e=chi2rnd(2*kE);%sampling residual variance
         si_a=ch_a/(2*lamA);si_p=ch_p/(2*lamD);si_e=ch_e/(2*lamE);
         va=(2*lamA)/ch_a;vp=(2*lamD)/ch_p;ve=(2*lamE)/ch_e;
    end
   if(iteration>=test) % entering the adaptive phase
        if(iteration==test)
                [co_mat]=co_va(fva,fvp,fve);%calculating the covariance matrix for the stored data
                
                               
        end
        iteration
        kA=1;kD=1;kE=1;
        lamA=0.001;lamD=0.001;lamE=0.001;
        
                old_v=1/2*log([va vp ve]);%values from tha current run
                new_v=old_v'+co_mat*randn(3,1);%sampling w*
                lg_jacb=-2*(new_v(1)-old_v(1)+new_v(2)-old_v(2)+new_v(3)-old_v(3));%Jacobian constant
                si_new=exp(-2*new_v);
                si_old=exp(-2*old_v);
                sn=mean(si_new);
                s_n=1/si_new(3);
		

                 %calculating likelihood ratio 

                new_cov=X*X'*1000000+A*(1/si_new(1))+D*(1/si_new(2))+eye(records)*(1/si_new(3));	%ZAZ'*sigma^2a+ZDZ'*sigma^2d+Isigma^2e=E for the new variance components 
                if(iteration==test)										 % here we have only one record per individual so ZAZ'=A and  eye function create an Identity matrix
                s_o=1/si_old(3);
                old_cov=X*X'*1000000+A*(1/si_old(1))+D*(1/si_old(2))+eye(records)*(1/si_old(3));
                old_cov=(1/s_o)*old_cov;
                old_cho=chol(old_cov,'lower');	
                old_det=2*sum(log(diag(old_cho)));
                sol_old=(old_cho\y);
                sol_old_t=sol_old'*sol_old;
                 end
		new_cov=(1/s_n)*new_cov;% for the current variance components
		new_cho=chol(new_cov,'lower');
		new_det=2*sum(log(diag(new_cho)));
        sol_new=(new_cho\y);
		sol_new_t=sol_new'*sol_new;
		new_lik_ratio=(records/2)*log(s_o)+0.5*old_det-(records/2)*log(s_n)-0.5*new_det-0.5*(1/s_n)*sol_new_t+0.5*(1/s_o)*sol_old_t;
                	
                  lg_mh_rat1=((kA-1)*log(si_new(1)))+(-lamA*si_new(1))+((kD-1)*log(si_new(2)))+(-lamD*si_new(2))+((kE-1)*log(si_new(3)))+(-lamE*si_new(3));
                lg_mh_rat2=((kA-1)*log(si_old(1)))+(-lamA*si_old(1))+((kD-1)*log(si_old(2)))+(-lamD*si_old(2))+((kE-1)*log(si_old(3)))+(-lamE*si_old(3));
              
               %lg_mh_rat1=-si_new(1)-si_new(2)-si_new(3);
                %  lg_mh_rat2=-si_old(1)-si_old(2)-si_old(3);
                lg_new_rat=(lg_mh_rat1-lg_mh_rat2)+new_lik_ratio +lg_jacb;		 %calculating the new ratio
                mh_rat=exp(lg_new_rat);
                   
                if(rand(1)< mh_rat)
                   va=1/si_new(1);
                   vp=1/si_new(2);
                   ve=1/si_new(3);
                   res=res+1;
                    s_o=1/si_new(3);
                    old_cov=X*X'*1000000+A*(1/si_new(1))+D*(1/si_new(2))+eye(records)*(1/si_new(3));
                    old_cov=(1/s_o)*old_cov;
                    old_cho=chol(old_cov,'lower');	
                    old_det=2*sum(log(diag(old_cho)));
                    sol_old=(old_cho\y);
                    sol_old_t=sol_old'*sol_old;
                  
                end
                                  
         
   end                              
     
            if(iteration<=test)
                    k_ratA=ve/va;   %sampling k ratio (of A)
                    k_ratE=ve/vp;   %sampling k ratio (of interaction factor)
   
                    if(k_ratA<1*10e-5)
                            k_ratA=1*10e-5;
                    end
                    if(k_ratE<1*10e-5)
                            k_ratE=1*10e-5;
                    end
                    
                    C = Cw + [bl11 bl12 bl13; bl21 A_inv*k_ratA bl23; bl31 bl32 D_inv*k_ratE];
            end

      %save every th round
      if(iteration > burnin)
         mn=mn+1;
         fva(mn) = va;
         fvp(mn) = vp;
         fve(mn) = ve;
      end
     
      if(iteration < burnin)
        bn=bn+1;
         fva_bn(bn) = va;
         fvp_bn(bn) = vp;
         fve_bn(bn) = ve;  
      end  
    end
   
cal=res/(iteration-test);
    
