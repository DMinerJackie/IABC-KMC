%%%%%%%%%%%%%%%%%%
%最大最小距离算法
%%%%%%%%%%%%%%%%%%
function addArray=testmaxmindistance915(x,ColonyNumber)
%这里的x相当于Conoly，addArray相当于Employed

% figure;
% for point=1:length(x)
%     plot(x(:,1),x(:,2),'r+');
% end
%centerNum=length(x)/2;%聚类中心点个数
centerNum=ColonyNumber;
index=0;%记录当前指针
[datarow,datacolumn]=size(x);%datarow表示数据的行数即样本数，datacolumn表示数据的属性个数



addArray=zeros(1,datacolumn);%新的聚类点加入的数组,初始化为一行，后面根据进入的聚类点数动态增加
maxArray=zeros(1,size(x,1));%存储最大值的数组
addArray(1,:)=x(1,:);%第一个聚类中心点赋给addArray数组的第一行


%***************寻找离第一个点最远的点****************
for k=1:size(x,1)
    %maxArray(1,k)=pdist([addArray(1,:);x(k,:)]');
    maxArray(1,k)=pdist2(addArray(1,:),x(k,:),'Euclidean');%9.20,这句话意思是说用欧氏距离算x和y之间的欧氏距离
end
[maxdistance,index]=max(maxArray);%maxdistance是求出的最大值，index是最大值所在的角标
    addArray(2,:)=x(index,:);%找到第二个聚类中心点并赋给第二行
    
    x(1,:)=[];
    maxArray(:,size(x,1))=[];
   % maxArray(:,index)=[];%对于本例数据，上句可以实现功能，但是如果最后一个是不是最大。。。没有必要了，这段代码只是找第二个
   % 聚类点，不需要循环效应。
    x(k-1,:)=[];%从x数组中去掉已作为中心点的点,因为前面有x(1,:)=[];本来是24行，现在是23行，本来找到的是18最大，现在要再减一，即k-1
    maxArray(:,size(x,1))=[];
    
 iformod=0;%计数器，用于后面mod(i,2)用，不取i是为了怕和后面的i混淆
    %*****************寻找第三到第centerNum个点******************
for m=1:centerNum-2%因为已经有两个点事先选好了，所以这里要在目标点数的基础上减去2
for i=1:size(x,1)
    [addrow,addcolumn]=size(addArray);
    [maxrow,maxcolumn]=size(maxArray);
    edArray=zeros(1,addrow);%用于存储x中每个元素到addArray中元素的距离
    %maxArray(length(x),:)=[];
for j=1:addrow
    %edArray(j)=pdist([addArray(j,:);x(i,:)]);
    edArray(j)=pdist2(addArray(j,:),x(i,:),'Euclidean');
end
maxArray(:,i)=max(edArray)*min(edArray);
end
iformod=iformod+1;
%if (mod(iformod,2)==1)%920为了防止陷入解全部在外围，所以采取这种方式规避
[maxforloop,index]=max(maxArray);
%else
 %   [maxforloop,index]=min(maxArray);
%end;
addArray((m+2),:)=x(index,:);
% disp('addArray的元素如下：');
% disp(addArray);
x(index,:)=[];
maxArray(:,size(x,1))=[];%每循环一次，x和maxArray的个数减少一个，因为需要将符合条件的点加入到addArray中，所以就少了一个
%另外，x和maxArray的关系是：x的行数与maxArray列数相等，可以利用maxArray中最大值的角标寻找x的行角标，从而知道该删哪一行
end
%     hold on
%     grid;
%     for point1=1:size(addArray,1)
%             plot(addArray(:,1),addArray(:,2),'b*');
%     end