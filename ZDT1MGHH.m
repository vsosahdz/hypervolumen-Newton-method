function [FH,DFH,DDFH,HH,DHH,DDHH]=ZDT1MGHH(n)
    x=sym('x',[n 1]);  
    assume(x,'real'); 
   
    f1=x(1);
    f2=(1.0+(9/(n-1))*(sum(x(2:end).^2)))*(2.0-sqrt(x(1)/(1.0+(9/(n-1))*(sum(x(2:end).^2))))); 
    
    h1=(x(1)-0.5)^2+(x(2)-0.4)^2-0.25;
    
    %Build the function and its jacobian and hessian matrices
    F=[f1;f2];
    DF=jacobian(F,x);
    DDF=[jacobian(DF(:,1),x);jacobian(DF(:,2),x)];
    
    H=h1;
    DH=jacobian(H,x);
    DDH=jacobian(DH,x);
    
    
    %Transform the function into matlab function to use them
    FH=matlabFunction(F,'vars',{x});
    DFH=matlabFunction(DF,'vars',{x}); 
    DDFH=matlabFunction(DDF,'vars',{x});
    
    HH=matlabFunction(H,'vars',{x});
    DHH=matlabFunction(DH,'vars',{x}); 
    DDHH=matlabFunction(DDH,'vars',{x});

end
