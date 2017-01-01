function ObjVal=Sphere(Colony,xd)
S=Colony.*Colony;%注意的Colony不是在runABC中的Colony，而是在调用ObjVal传过来的Employed
ObjVal=sum(S');

