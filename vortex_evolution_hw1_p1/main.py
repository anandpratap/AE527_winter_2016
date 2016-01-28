from pylab import *
import fortutils as fu
import time
class VortexSolver(object):
    def __init__(self):
        self.nx = 150
        self.ny = 75
        self.L = 10.0
        self.xr = np.linspace(-self.L, self.L, self.nx)
        self.yr = np.linspace(-self.L/2, self.L/2, self.ny)
        self.x, self.y = np.meshgrid(self.xr, self.yr)
        self.x = self.x.T
        self.y = self.y.T
        self.alpha = np.zeros([self.nx, self.ny])
        self.omega = np.zeros([self.nx, self.ny])
        self.u = np.zeros([2, self.nx, self.ny])
        self.init_omega()
        self.init_strength()
        self.write_freq = 10
        self.realtimeplot = False
        self.fortran = True
        self.integrator = 4
        
    def init_omega(self):
        Tau = 1.0
        rc = 1.0
        
        xc, yc = -5.0, 0.0
        r = np.sqrt((self.x - xc)**2 + (self.y - yc)**2)
        self.omega += Tau/pi/rc**2*(exp(-r**2/rc**2))

        xc, yc = 5.0, 0.0
        r = np.sqrt((self.x - xc)**2 + (self.y - yc)**2)
        self.omega += -Tau/pi/rc**2*(exp(-r**2/rc**2))

        Tau = 0.1
        rc = 0.5
        xc, yc = -3.0, 0.0
        r = np.sqrt((self.x - xc)**2 + (self.y - yc)**2)
        self.omega -= Tau/pi/rc**2*(exp(-r**2/rc**2))

        xc, yc = 3.0, 0.0
        r = np.sqrt((self.x - xc)**2 + (self.y - yc)**2)
        self.omega += Tau/pi/rc**2*(exp(-r**2/rc**2))
        
    def init_strength(self):
        dx = self.xr[1] - self.xr[0]
        dy = self.yr[1] - self.yr[0]
        A = dx*dy
        self.alpha = self.omega*A

    def calc_velocity(self):
        delta = 1
        for i in range(self.nx):
            for j in range(self.ny):
                dx = self.x[i,j] - self.x
                dy = self.y[i,j] - self.y
                r2 = dx**2 + dy**2
                rt2 = r2/delta**2
                ux_i = dy*self.alpha/(2.0*pi*(r2 + 1e-16))*(1.0 - exp(-rt2))
                uy_i = -dx*self.alpha/(2.0*pi*(r2 + 1e-16))*(1.0 - exp(-rt2))
                self.u[0,i,j] = np.sum(ux_i)
                self.u[1,i,j] = np.sum(uy_i)

    def write_tecplot(self, step):
        step = str(step)
        namestr = (5-len(step))*'0' + step
        f = open("data/step_%s.tec"%namestr, "w")
        f.write("""Variables="X","Y","OMEGA","ALPHA", "U", "V"\nZone I=     %i,J=     %i, F=POINT\n"""%(self.nx, self.ny))
        for i in range(self.nx):
            for j in range(self.ny):
                f.write("%f %f %f %f %f %f\n"%(self.x[i,j], self.y[i,j], self.omega[i,j], self.alpha[i,j], self.u[0,i,j], self.u[1,i,j]))
        f.close()
        

    def run(self):
        t = 0.0
        step = 0
        tf = 1000.0
        dt = 1e-1
        ion()

        while t < tf:
            start = time.time()
            if self.fortran:
                if self.integrator == 1:
                    self.u = fu.calc_velocity(self.x, self.y, self.alpha)
                    self.x += self.u[0,:,:]*dt
                    self.y += self.u[1,:,:]*dt

                elif self.integrator == 4:
                    x = self.x.copy()
                    y = self.y.copy()
                    k1 = fu.calc_velocity(x, y, self.alpha)
                    x = self.x + k1[0,:,:]*dt/2.0
                    y = self.y + k1[1,:,:]*dt/2.0

                    k2 = fu.calc_velocity(x, y, self.alpha)
                    x = self.x + k2[0,:,:]*dt/2.0
                    y = self.y + k2[1,:,:]*dt/2.0
                    
                    k3 = fu.calc_velocity(x, y, self.alpha)
                    x = self.x + k3[0,:,:]*dt
                    y = self.y + k3[1,:,:]*dt

                    k4 = fu.calc_velocity(x, y, self.alpha)
                    self.u[0,:,:] = (k1[0,:,:] + 2.0*k2[0,:,:] + 2.0*k3[0,:,:] + k4[0,:,:])/6.0
                    self.u[1,:,:] = (k1[1,:,:] + 2.0*k2[1,:,:] + 2.0*k3[1,:,:] + k4[1,:,:])/6.0
                    self.x += self.u[0,:,:]*dt
                    self.y += self.u[1,:,:]*dt
                else:
                    raise ValueError("Wrong integrator.")
            else:
                self.calc_velocity()
                self.x += self.u[0,:,:]*dt
                self.y += self.u[1,:,:]*dt

            cpu_dt = time.time() - start
            
            if step%self.write_freq == 0:
                self.write_tecplot(step)
                if self.realtimeplot:
                    figure(1)
                    clf()
                    contourf(self.x, self.y, self.alpha, 50)
                    colorbar()
                    figure(2)
                    clf()
                    plot(self.x, self.y, 'ro')
                    pause(0.0001)

            t += dt
            step += 1
            print "TIME: %.8e, CPU TIME: %.8fs"%(t,cpu_dt)
                            
        
if __name__ == "__main__":
    vs = VortexSolver()
    vs.realtimeplot = True
    vs.write_freq = 10
    vs.run()
