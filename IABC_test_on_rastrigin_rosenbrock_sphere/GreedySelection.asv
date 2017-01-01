
function [Colony Obj Fit oBas]=GreedySelection(Colony1,Colony2,ObjEmp,ObjEmp2,FitEmp,FitEmp2,fbas,ABCOpts,i)

oBas=fbas;
Obj=ObjEmp;
Fit=FitEmp;
Colony=Colony1;
if (nargin==8)%Inside the body of a user-defined function, NARGIN returns the number of input arguments that were used to call the function. 
for ind=1:size(Colony1,1)%ind=1：5，对所有食物源进行贪婪选择
    if (FitEmp2(ind)>FitEmp(ind))%如果Vi的适应度值大于Xi的，替换，，
        oBas(ind)=0;
         %zj因为这是已经被新的位置更新了，所以其开采度应该置为零，表示这是第一次，没有被开采过
        Fit(ind)=FitEmp2(ind);
        Obj(ind)=ObjEmp2(ind);
        Colony(ind,:)=Colony2(ind,:);
    else%否则不变，并且计数器bas+1
        oBas(ind)=fbas(ind)+1;
         %zj因为新的位置的适应度没有当前的好（大），所以在当前位置上仍保留当前解，表示当前又被开采了一次
        Fit(ind)=FitEmp(ind);
        Obj(ind)=ObjEmp(ind);
        Colony(ind,:)=Colony1(ind,:);
    end;
end; %for
end; %if
if(nargin==9)%第i个引领蜂被跟随，只对第i个食物源进行贪婪选择
    ind=i;
    if (FitEmp2(ind)>FitEmp(ind))
        oBas(ind)=0;
        Fit(ind)=FitEmp2(ind);
        Obj(ind)=ObjEmp2(ind);
        Colony(ind,:)=Colony2(ind,:);
    else
        oBas(ind)=fbas(ind)+1;
        Fit(ind)=FitEmp(ind);
        Obj(ind)=ObjEmp(ind);
        Colony(ind,:)=Colony1(ind,:);
    end;
end; 
    