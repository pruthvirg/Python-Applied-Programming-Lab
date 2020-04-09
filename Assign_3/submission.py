from pylab import *
from scipy.special import *
import numpy as np

f = open('./fitting.dat','r')
rows = f.readlines()
no_backslashn = []
for i in rows:
    no_backslashn.append(i.split('\n')[0])

matrix = []
for i in no_backslashn:
    matrix.append(i.split(' '))

#extracting the signals
signal_dict = {'sig0':[],
'sig1' : [],
'sig2' : [],
'sig3' : [],
'sig4' : [],
'sig5' : [],
'sig6' : [],
'sig7' : [],
'sig8' : [],
'sig9' : []}

for i in range(0,10):
    for t in matrix:
        key = 'sig' + str(i)
        signal_dict[key].append(float(t[i]))

def noise(sig):
    return np.std(np.asarray([x1 - x2 for (x1, x2) in zip(sig,signal_dict['sig9'])]))

noise_list = []
for i in range(1,10):
    sig = "sig" + str(i)
    noise_list.append(noise(signal_dict[sig]))

M = np.matrix(matrix)
plot(M[:,0],M[:,1:9])
value=('True value','σ = {:.3f}'.format(noise_list[0]),'σ = {:.3f}'.format(noise_list[1]),'σ = {:.3f}'.format(noise_list[2]),'σ = {:.3f}'.format(noise_list[3]),'σ = {:.3f}'.format(noise_list[4]),'σ = {:.3f}'.format(noise_list[5]),'σ = {:.3f}'.format(noise_list[6]),'σ = {:.3f}'.format(noise_list[7]),'σ = {:.3f}'.format(noise_list[8]))
gca().legend(value)
xlabel(r'$t$',size=20)
ylabel(r'$f(t)+noise$',size=20)
title(r'Data to be fit')
grid(True)
show()

def err(l):
    err = 0
    for i in range(len(l)-1):
        err = err + (l[i] - signal_dict['sig9'][i])**2
    return err

for i in range(1,10):
    sig_name = "sig" + str(i)
    signal_dict[i]= err(signal_dict[sig_name]) 


#code for plotting bessel functions with given params

def plot_bessels(a,b):
    x_points = np.asarray([float(i) for i in signal_dict["sig0"]])
    y_bessel = float(a)*jv(0,x_points)+float(b)*x_points
    title(r'True function')
    xlabel(r'$t$',size=20)
    ylabel(r'$f(t)$',size=20)
    plot(x_points,y_bessel)

plt.plot(linspace(0,10,101),signal_dict["sig9"])
plt.grid(True)
errorbar(linspace(0,10,21),signal_dict['sig1'][::5],0.096,fmt='ro')
gca().legend(('True value','error bar'))
xlabel(r'$t$',size=20)
title(r'Data points with σ=0.096 along with exact function')
grid(True)
show()

#problem 6:
def g2(t,a,b):
    y = jv(2,t)
    #print(y)
    m = c_[y,t]
    a_b  = np.transpose(np.asarray([a,b]))
    return np.matmul(m,a_b)

'''
sol for question 7:
'''

A = linspace(0,2,21)
B = linspace(-0.2,0,21)

f_true = g2(signal_dict["sig0"],1.05,-0.105)

def cal_error(a,b):
    return(np.mean(np.square(f_true-g2(signal_dict["sig0"],a,b))))

err_matrix = np.zeros((len(A),len(A)))
for i in range(len(A)):
    for t in range(len(B)):
        err_matrix[i,t] = cal_error(A[i],B[t])

#plotting the contour of the epsilon matrix
CS=contour(A,B,err_matrix,10)
plot(1.05,-0.105,'ro',)
title(r'Contour plot of epsilon_ij')
matplotlib.pyplot.text(1.05, -0.105, "Exact location")
xlabel(r'$A$',size=20)
ylabel(r'$B$',size=20)
clabel(CS, inline=1, fontsize=10)
show()

#finding optimal value by linalg.lstsq function
def find_optimum_a_b(signal):
    y = jv(2,signal)
    m = c_[y,signal]
    m_true = f_true
    a, b = np.linalg.lstsq(m, f_true)[0]
    return [a,b]

find_optimum_a_b(signal_dict["sig1"])

optimal_a = []
optimal_b = []
for i in range(1,10):
    sig = "sig" + str(i)
    optimal_a.append(np.square(find_optimum_a_b(signal_dict[sig])[0]-1.05)/40)
    optimal_b.append(np.square(find_optimum_a_b(signal_dict[sig])[1]+0.105)/40)

#plotting error in linear scale
plot(list(reversed(noise_list)),list(reversed(optimal_a)) , 'ro',linestyle='--',dashes=(2, 10))
plot(list(reversed(noise_list)),list(reversed(optimal_b)), 'go',linestyle='--',dashes=(2, 10))
gca().legend(('Aerr','berr'))
xlabel(r'$Noise stdev$',size=15)
ylabel(r'$Ms Error$',size=15)
title(r'Plot of error in linear scale',size=18)
grid(True)
show()

#plotting error in loglog scale
loglog(list(reversed(noise_list)),list(reversed(optimal_a)) , 'ro')
loglog(list(reversed(noise_list)),list(reversed(optimal_b)) , 'go')
gca().legend(('Aerr','berr'))
xlabel(r'$Noise stdev$',size=15)
ylabel(r'$Ms Error$',size=15)
title(r'variation of error in loglog scale',size=15)
grid(True)
show()
    
