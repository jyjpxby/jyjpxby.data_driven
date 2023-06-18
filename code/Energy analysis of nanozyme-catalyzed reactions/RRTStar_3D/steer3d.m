function A = steer3d(qr, qn, val, eps)
data1=importdata('steer3d_Settings.txt');%Each time you make a change, replace the value in this file with the corresponding value in Table S4.
data2=data1.data;
x_min = data2(1,1);
y_min = data2(1,2);
z_min = data2(1,3);
   qnew = [data2(1,12)-x_min data2(1,13)-y_min];
   if val >= eps
       qnew(1) = qn(1) + ((qr(1)-qn(1))*eps)/dist_3d(qr,qn);
       qnew(2) = qn(2) + ((qr(2)-qn(2))*eps)/dist_3d(qr,qn);
       qnew(3) = qn(3) + ((qr(3)-qn(3))*eps)/dist_3d(qr,qn);
   else
       qnew(1) = qr(1);
       qnew(2) = qr(2);
       qnew(3) = qr(3);
   end
   
   A = [qnew(1), qnew(2), qnew(3)];
end