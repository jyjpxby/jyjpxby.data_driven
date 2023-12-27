function A = steer3d(qr, qn, val, EPS)
x_min = -15;
y_min = -15;
   qnew = [0-x_min 0-y_min];
   if val >= EPS
       qnew(1) = qn(1) + ((qr(1)-qn(1))*EPS)/dist_3d(qr,qn);
       qnew(2) = qn(2) + ((qr(2)-qn(2))*EPS)/dist_3d(qr,qn);
       qnew(3) = qn(3) + ((qr(3)-qn(3))*EPS)/dist_3d(qr,qn);
   else
       qnew(1) = qr(1);
       qnew(2) = qr(2);
       qnew(3) = qr(3);
   end
   
   A = [qnew(1), qnew(2), qnew(3)];
end