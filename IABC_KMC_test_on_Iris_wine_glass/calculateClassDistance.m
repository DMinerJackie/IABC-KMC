function [newCenter,class,classCounterDistance]=calculateClassDistance(x,center,datarow,centerNum,sum1,kindNum)

if(nargin==4)
  for i=1:datarow
    for j=1:centerNum
      distancedata(j,i)=pdist2(x(:,i)',center(:,j)','Euclidean');%920,明确采用的是欧氏距离计算两个点之间的距离
        % distancedata(j,i)=pdist([x(:,i),center(:,j)]');
  %   distancedata(j,i)=(x(:,i)-center(:,j))'*(x(:,i)-center(:,j));%计算各第i点到就各个中心的距离
    end
  temp=find(distancedata(:,i)==min(distancedata(:,i)));%9.17 之所以这么写是因为怕出现如x=[1,1,3],find(x==min(x))会返回两个值而不符合实际情况，
  %所以现将其以数组形式存储起来，再用下面的t=temp(1,1)取第一个就可以了
    t=temp(1,1);%i到第t个中心近
    %**************9.17**************
    %t就是这个点到所有聚类中心中距离最短的那个中心点所代表的的类，这里有三个中心点，如果t==2，表示i这个点到第二个中心点的距离最短，所以这个
    %点也应该归属到第二类中，即t==2
    class(i)=t;
    classCounterDistance(i,:)=[min(distancedata(:,i)),t];%这个将这个i点对应从属于哪个类以及到这个类的聚类存储起来，
    %这个是为了以后计算准则函数用的，就是把同类的距离都加起来
    %**************9.17**************
  end
  newCenter=center;
end

if(nargin==6)
  for i=1:datarow
    for j=1:centerNum
      distancedata(j,i)=pdist2(x(:,i)',center(:,j)','Euclidean');%920,明确采用的是欧氏距离计算两个点之间的距离
        % distancedata(j,i)=pdist([x(:,i),center(:,j)]');
  %   distancedata(j,i)=(x(:,i)-center(:,j))'*(x(:,i)-center(:,j));%计算各第i点到就各个中心的距离
    end
  temp=find(distancedata(:,i)==min(distancedata(:,i)));%9.17 之所以这么写是因为怕出现如x=[1,1,3],find(x==min(x))会返回两个值而不符合实际情况，
  %所以现将其以数组形式存储起来，再用下面的t=temp(1,1)取第一个就可以了
    t=temp(1,1);%i到第t个中心近
    %**************9.17**************
    %t就是这个点到所有聚类中心中距离最短的那个中心点所代表的的类，这里有三个中心点，如果t==2，表示i这个点到第二个中心点的距离最短，所以这个
    %点也应该归属到第二类中，即t==2
    class(i)=t;
    classCounterDistance(i,:)=[min(distancedata(:,i)),t];%这个将这个i点对应从属于哪个类以及到这个类的聚类存储起来，
    %这个是为了以后计算准则函数用的，就是把同类的距离都加起来
    %**************9.17**************
    
%     for j=1:length(x)/2
%         if t==j
%             sum1(:,j)=sum1(:,j)+x(:,i);
%             kindNum(j)=kindNum(j)+1;
%         end
%     end
sum1(:,t)=sum1(:,t)+x(:,i);
kindNum(t)=kindNum(t)+1;
%这样写可以省去不必要的循环，显然t和j可以一一对应，直接用t替换j

  end
   
  for t=1:centerNum%924将j全部替换成t
    newCenter(:,t)=sum1(:,t)./kindNum(t);%计算新的中心
 end
%920
end