# Question 2
u=[2,-4,0]
v=[3,1.5,-7]

# 2a)
w=2*u + 5*v
# 2b) dot product
d=dot(u,v)
# 2c) magnitude of u is the norm()
l=norm(u)
# 2d) divide u by the magnitude to get the unit vector u1
u1=u ./ norm(u)
# 2e) cross product will generate a vector orthogonal to u and v
n=cross(u,v)
# 2f) same as 2c the dot product shows the projection of u along v
proj_u_along_v = dot(u,v)

# Question 3
A=[1,-3,5;2,-4,3;0,1,-1]
B=[1,-1,0,0;-3,0,7,-6;2,1,-2,-1]
I3=eye(3)

# 3a) determinant = 5
d=det(A)
# 3b)
C=2*A + 4*I3
# 3c)
D=inv(A)
# 3d) can invert B, a non-square matrix can't be inverted
# E=inv(B)
# 3e) B*A can't be multiplied, if B is size pxq then A must be qxr
# B is 3x4,  A is 3x3
# F=B*A
# 3f)
G=(A*B)'
# 3g) B transponsed is 4x3, A transposed is 3x3, result is 4x3
H=B'*A'
### (A*B)' = B'*A'
[narz@centos7 octave_projects]$ cat lab1.m
# Question 2
u=[2,-4,0]
v=[3,1.5,-7]

# 2a)
w=2*u + 5*v
# 2b) dot product
d=dot(u,v)
# 2c) magnitude of u is the norm()
l=norm(u)
# 2d) divide u by the magnitude to get the unit vector u1
u1=u ./ norm(u)
# 2e) cross product will generate a vector orthogonal to u and v
n=cross(u,v)
# 2f) same as 2c the dot product shows the projection of u along v
proj_u_along_v = dot(u,v)

# Question 3
A=[1,-3,5;2,-4,3;0,1,-1]
B=[1,-1,0,0;-3,0,7,-6;2,1,-2,-1]
I3=eye(3)

# 3a) determinant = 5
d=det(A)
# 3b)
C=2*A + 4*I3
# 3c)
D=inv(A)
# 3d) can invert B, a non-square matrix can't be inverted
# E=inv(B)
# 3e) B*A can't be multiplied, if B is size pxq then A must be qxr
# B is 3x4,  A is 3x3
# F=B*A
# 3f)
G=(A*B)'
# 3g) B transponsed is 4x3, A transposed is 3x3, result is 4x3
H=B'*A'
### (A*B)' = B'*A'

close all;
x=linspace(-2*pi, 2*pi, 1000);
y=x.^2 .* sin(x);
figure(1); 
plot(x,y, 'r-', 'linewidth', 2); hold;

y2=x.^2;
y3=-1.*x.^2;
plot(x, y2, 'k:', x, y3, 'k:');
legend("x^2*sin(x)", "x^2", "-x^2");
grid on;
title('Lab1');
