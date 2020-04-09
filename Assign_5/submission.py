from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3


Nx = 25        
Ny = 25         
radius = 0.35
Niter = 1500   
errors = np.zeros(Niter)


x = np.linspace(-0.5,0.5,25)    
y = np.linspace(0.5,-0.5,25)   
X,Y = meshgrid(x,y)            
phi = np.zeros((Nx,Ny))         
ii = where(X*X + Y*Y <= radius*radius) 
phi[ii] = 1.0                  

contour(X,Y,phi)
plot(x[ii[0]],y[ii[1]],'ro')
grid()
title('Contour plot of initial potential')
xlabel('x')
ylabel('y')
show()

newphi = np.zeros((Nx,Ny)) 
for k in range(Niter):
    oldphi = phi.copy()   
    newphi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2] + phi[1:-1,2:] + phi[0:-2,1:-1] + phi[2:,1:-1])  
    
    newphi[1:-1,0] = newphi[1:-1,1]       
    newphi[1:-1,Nx-1] = newphi[1:-1,Nx-2]
    newphi[0,1:-1] = newphi[1,1:-1]
    newphi[ii] = 1.0
    
    errors[k] = max(np.absolute(np.subtract(oldphi.flatten(),newphi.flatten())))   
    phi = newphi.copy() 



xError = np.linspace(1,Niter,1500)  
yError = np.log(errors)              
A=np.zeros((Niter,2))               
A[:,0] = 1
A[:,1] = xError
const = lstsq(A,yError)[0]           
yError = const[0] + const[1]*xError  
yError = np.exp(yError)

semilogy(xError,errors)
show()

loglog(np.arange(1,1501,50),errors[0::50],'ro')
loglog(xError,errors)
show()

xError2 = np.linspace(501,Niter,1000)
yError2 = np.log(errors[500:])
B=np.zeros((Niter-500,2))
B[:,0] = 1
B[:,1] = xError2
const = lstsq(B,yError2)[0]
yError2 = const[0] + const[1]*xError2
yError2 = np.exp(yError2)


semilogy(np.arange(1,1501,50),errors[0::50],'ro')
plot(xError,yError)
plot(xError2, yError2)
grid()
title('Error plot')
xlabel('No. of iterations')
ylabel('Error')
legend(('Calculated Error','Fit 1 (all iterations)','Fit 2 (>500 iterations)'))
show()




fig1 = figure(4)
ax = p3.Axes3D(fig1)
title('The 3-D surface plot of the potential')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Potential $(\phi)$')
surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
show()


contour(x,y,phi)
plot(x[ii[0]],y[ii[1]],'ro')
xlabel('x')
ylabel('y')
title('Contour plot of final potential')
grid()
show()



Jx = np.zeros((Nx,Ny))
Jy = np.zeros((Nx,Ny))

Jy[1:-1,1:-1] = 0.5*(phi[1:-1,2:] - phi[1:-1,0:-2])
Jx[1:-1,1:-1] = 0.5*(phi[2:,1:-1] - phi[0:-2,1:-1])




plot(x[ii[0]],y[ii[1]],'ro')
xlabel('x')
ylabel('y')
title('Vector plot of the current flow')
quiver(y,x,Jy[::-1,:],Jx[::-1,:])
contour(x,y,phi)
show()

