way=input('automatic generation of initial points 1 or manual 2?');
cv=input('Trace to leave out? 0 if none, number of the trace to leave out or vector if several traces to leave out');
area=input('Hippocampus 1 or neocortex 2?');
if area==1
    set=input('Subthreshold protocol 1 or STDP protocol 2?');
elseif area==2
    set=input('Letzkus 1 or Sjostrom 2?');
end

init=xlsread('valeurs_init.xlsx',1); %excel file, each row contains one set of initial parameters
lb=[2;2;5;2.5;10^-5;10^-5;2;50;2];
ub=[30;60;30;15;10^-2;10^-2;60;50000;100];

if way==1
    
    for c=1:length(init)
        
        A=[0 0 -1 1 0 0 0 0 0];
        b=0;
        Aeq=[0 0 0 0 0 0 0 0 0];
        beq=0;
        y0=init(c,:);
        
        if area==1
            [y,STDPval]=fmincon(@(y)model_hipp(y,cv,set),y0,A,b,Aeq,beq,lb,ub);
        elseif area==2
            [y,STDPval]=fmincon(@(y)model_neocortex1(y,set),y0,A,b,Aeq,beq,lb,ub);
        end
        y_store(c,:)=y;
        error_store(c)=STDPval;
        
        clear y STDPval
        
    end
    save(strcat('paramfit_area',num2str(area),'exp',num2str(set),num2str(cv)),'y_store','error_store','lb','ub')
    
elseif way==2
    rng default % For reproducibility
    gs = GlobalSearch;
    if area==1
        lse=@(y)model_hipp(y,cv,set);
    elseif area==2
        lse = @(y)model_neocortex1(y,set);
    end
    for it=1:length(init)
        xval=init(it,:);
        
        problem = createOptimProblem('fmincon','x0',xval,'objective',lse,'lb',lb,'ub',ub,'Aineq',[0 0 -1 1 0 0 0 0 0],'bineq',0);
        [x,fval,exitflag,output,solutions] = run(gs,problem);
        
        Final_out{it}.x=x;
        Final_out{it}.fval=fval;
        Final_out{it}.exitflag=exitflag;
        Final_out{it}.output=output;
        Final_out{it}.solutions=solutions;
    end
    save(strcat('paramfit_area',num2str(area),'exp',num2str(set),num2str(cv)),'Final_out')
end

