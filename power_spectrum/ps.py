# create some "fake" 3-d multimode data and create a power spectrum to
# illustrate that we have the scaling and normalization correct

import numpy as np
import pylab

N = 64

index = -2.0

xmin = ymin = zmin = 0.0
xmax = ymax = zmax = 10.0

L = np.array([xmax - xmin, ymax - ymin, zmax - zmin])

phi = np.zeros((N,N,N), dtype=np.float64)

x = (np.arange(N)+0.5)*(xmax - xmin)/N + xmin
y = (np.arange(N)+0.5)*(ymax - ymin)/N + ymin
z = (np.arange(N)+0.5)*(zmax - zmin)/N + zmin

x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

modes = 16

A_0 = 1000.0


# first pass -- count how many hits each k gets
weights = {}

for m in range(1,modes+1):
    for n in range(1,modes+1):
        for p in range(1,modes+1):

            ii2 = m**2 + n**2 + p**2
            
            if not ii2 in weights.keys():
                weights[ii2] = 1
            else:
                weights[ii2] += 1


# amplitude of smallest k
print "min{k_x} = ", 1/L[0]
print "amplitude of smallest k = ", A_0*np.sqrt(1/L[0]**2 + 1/L[1]**2 + 1/L[2]**2)**index/weights[3]


# compute the function we will find the power spectrum of
for m in range(1,modes+1):
    k_x = m/L[0]
    print "m = ", m

    for n in range(1,modes+1):
        k_y = n/L[1]

        for p in range(1,modes+1):
            k_z = p/L[2]

            ii2 = m**2 + n**2 + p**2

            k = np.sqrt(k_x**2 + k_y**2 + k_z**2)
            A = A_0*k**index/weights[ii2]
            
            phi += A*np.sin(k_x*x3d + k_y*y3d + k_z*z3d)


# now do the power spectrum
phi_hat = np.fft.fftn(phi) #[0:N/2+1,0:N/2+1,0:N/2+1]
#phi_hat = 8.0*phi_hat/N**3
phi_hat = phi_hat/N**3

phi_hat = abs(phi_hat)**2

# Parseval's theorem: sum of |phi(x)|**2 = sum of |phi_hat(k)|**2
print "sum of |phi(x)|**2    = {}".format(np.sum(phi))
print "sum of |phihat(k)|**2 = {}".format(np.sum(phi_hat))
print "DC Offset = {}".format(phi_hat[0,0,0])

kx = np.fft.fftfreq(N)[0:N/2+1]
kx[-1] *= -1
kx = kx*N/L[0]

ky = np.fft.fftfreq(N)[0:N/2+1]
ky[-1] *= -1
ky = ky*N/L[2]

kz = np.fft.fftfreq(N)[0:N/2+1]
kz[-1] *= -1
kz = kz*N/L[2]

kmin = 0
kmax = np.sqrt(np.max(kx)**2 + np.max(ky)**2 + np.max(kz)**2)

num = np.floor(np.sqrt(N**2 + N**2 + N**2))

bins = np.linspace(kmin, kmax, num+1)

kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")

k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

whichbin = np.digitize(k.flat, bins)

ncount = np.bincount(whichbin)

E_spectrum = np.zeros(len(ncount)-1, dtype=np.float64)

for n in range(ncount):
    if not ncount[n] == 0: 
        E_spectrum[n-1] = np.sum(phi_hat.flat[whichbin==n]) #/ncount[n]


k = bins[1:n]
E_spectrum = E_spectrum[0:len(k)]

pylab.loglog(k, E_spectrum)

ii = np.argmax(E_spectrum)
kmax = k[ii]
Emax = E_spectrum[ii]

pylab.loglog(k, Emax*(k/kmax)**index, ls=":", color="0.5")

pylab.savefig("ps.png")


