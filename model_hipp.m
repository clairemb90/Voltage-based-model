function S=model_hipp(y,cv,set)

dt=0.1;

int=1000; %interval %ms
inter=length([dt:dt:1000]);

[trace,filetoread,raw]=xlsread('traces_hipp',set);

if set==1
    supra=[20,0,0,0,19,1,0,16,0,3,5,0];
    nb_pair=1;
    data=15;
    com=[1 2 3 4 5 6 7 8 9 10 11 12];
    trans=[1 3 4 5 6 8 9 10 12 13 14 15];
elseif set==2
    data=4;
    nb_pair=50;
    com=[1 2 3 4];
    trans=[1 3 4 5];
    even(1,:)=[21 1 31 36 50 22 23 10 13 19 6 8 45 40 29 43 14 46 47 25 49 32 5 17 11 3 37 35 30 44 26 7 33 42 18 20 12 16 34 2 38 4 28 15 24 41 48 27 39 9];
end


dw=zeros(data,1);
dw_ltd=zeros(data,1);
dw_ltp=zeros(data,1);
nb_tot=nb_pair*int;

tspan=[dt:dt:nb_tot]; % ms
tproto=length(tspan);

Ed=0;
theta_moins_ltp=y(3);
theta_moins_ltd=y(4);
ce_sup=[];

for bou=1:length(cv)
    tenpo=com-cv(bou);
    ce_sup=find(tenpo~=0);
    com=com(ce_sup);
    trans=trans(ce_sup);
end

ww=0.5*ones(data,1);

for expdata=1:data
    u=zeros(tproto,1);
    titi=zeros(tproto,1);
    toto=zeros(tproto,1);
    wi=zeros(tproto,1); %weight
    if set==1
        for pp=0:nb_pair-1
            if nb_pair==1
                u(50+pp*inter:50+pp*inter+length(trace)-12)=trace(12:length(trace),expdata)-mean(trace(1:5,expdata));
                u(isnan(u)) = 0;
                if expdata==14 || expdata==15
                    titi(pp*inter+50+40/dt)=1;
                else
                    titi(pp*inter+50)=1;
                end
            else
                if  com(expdata)==1 || com(expdata)==5 || com(expdata)==8
                    if even(com(expdata),pp+1)>supra(com(expdata))
                        u(50+pp*inter:50+pp*inter+length(trace)-12)=trace(12:length(trace),trans(expdata))-mean(trace(1:5,trans(expdata)));
                        u(isnan(u)) = 0;
                    else
                        u(50+pp*inter:50+pp*inter+length(trace)-12)=trace(12:length(trace),trans(expdata)+1)-mean(trace(1:5,trans(expdata)+1));
                        u(isnan(u)) = 0;
                    end
                else
                    u(50+pp*inter:50+pp*inter+length(trace)-12)=trace(12:length(trace),trans(expdata))-mean(trace(1:5,trans(expdata)));
                    u(isnan(u)) = 0;
                end
                if com(expdata)==11 || com(expdata)==12
                    titi(pp*inter+50+40/dt)=1;
                else
                    titi(pp*inter+50)=1;
                end
            end
        end
    elseif set==2
        for pp=0:nb_pair-1
            if com(expdata)==1
                if even(com(expdata),pp+1)>20
                    u(50+pp*inter:50+pp*inter+length(trace)-1)=trace(:,trans(expdata));
                    u(isnan(u)) = 0;
                else
                    u(50+pp*inter:50+pp*inter+length(trace)-1)=trace(:,trans(expdata)+1);
                    u(isnan(u)) = 0;
                end
            else
                u(50+pp*inter:50+pp*inter+length(trace)-1)=trace(:,trans(expdata));
                u(isnan(u)) = 0;
            end
            titi(pp*inter+50)=1;
        end
    end
    
    x_barre=zeros(tproto,1);%a presynaptic spike increases the trace x of some biophysical quantity, then decays exp.
    ltd_veto=zeros(tproto,1);
    ltd_veto(1)=0;
    u_bare=zeros(tproto,1); %low-pass filtered version of the postsynaptic membrane potential
    u_bare_ltd=zeros(tproto,1);
    
    x_barre(1)=0;
    u_bare(1)=Ed;
    u_bare_ltd(1)=Ed;
    
    wi(1)=0.5;   
    tsim=2;
    
    timeexp=[];
    while tsim<=tproto
        
        if titi(tsim)==1
            timeexp=[timeexp tsim];
        end
        sum_exp=0;
        for kti=1:length(timeexp)
            sum_exp=sum_exp+exp((timeexp(kti)-tsim)*dt/y(1));
        end
        
        x_barre(tsim)=sum_exp;
        
        u_bare(tsim)=u_bare(tsim-1)+dt*(-u_bare(tsim-1)+u(tsim-1))/y(2);
        u_bare_ltd(tsim)=u_bare_ltd(tsim-1)+dt*(-u_bare_ltd(tsim-1)+u(tsim-1))/y(7);
        
        AA=u_bare(tsim-1)>theta_moins_ltp;
        BB=u_bare_ltd(tsim-1)>theta_moins_ltd(tsim-1);
        
        toto(tsim-1)=y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp); 
        ltd_veto(tsim)=ltd_veto(tsim-1)+dt*(-ltd_veto(tsim-1)+y(8)*(toto(tsim-1)))/y(9);
        theta_moins_ltd(tsim)=y(4)+ltd_veto(tsim);
        
        %final rule, linear
        wi(tsim)=wi(tsim-1)+y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp)*dt-y(6)*x_barre(tsim-1)*BB*(u_bare_ltd(tsim-1)-theta_moins_ltd(tsim-1))*dt;
                
        tsim=tsim+1;
        
    end
    
    dw(expdata)=wi(length(wi));
    
end


if set==1
    for expdata=1:length(trans)
        if trans(expdata)==1 || trans(expdata)==6 || trans(expdata)==10
            norm_weight(expdata)=1+(60-supra(com(expdata)))*(dw(trans(expdata))-ww(1))/0.5 + supra(com(expdata))* (dw(trans(expdata)+1)-ww(1))/0.5;
        else
            norm_weight(expdata)=(0.5+60.*(dw(trans(expdata))-ww(1)))./ww(1);
        end
    end
    norm_weight=norm_weight';
elseif set==2
    norm_weight=(dw./ww);
end

if set==1
    plas=[122;101.4;95.5;100;131;96.6;100;119.3;104.5;104.3;60;60]/100;
    sd_all=[1;1;1;1;1;1;1;1;1;1;2;2];
elseif set==2
    plas=[1.45; 1.04;1.02;1.02];
    sd_all=[1; 1;1;1];
end
norm_w_exp=[plas(com,1)];
sd=sd_all(com,1);

S=sum((norm_weight-norm_w_exp).^2./sd.^2)

end