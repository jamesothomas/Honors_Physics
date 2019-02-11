# this is where the numerical integration is done
# Be careful changing any of this

class slinky(object):
    
    def __init__(self, NX=40, L=0.5, K=5.0, M=1.0, g=9.8, NT=1000, T = 1.0, held=False):
        self.NX = NX
        self.NT = NT
        self.T = T
        self.dt = self.T / self.NT
        self.L = L 
        self.K = K
        self.M = M
        self.G = g
        self.held = held
        self.sk = self.K * (self.NX - 1)
        self.sl = self.L / (self.NX - 1)
        self.mm = self.M / self.NX

    def get_xdot(self, x):
        '''
        Note for whoever tries to interpret this later, It uses 1 list "xdot" to hold both xdot as well as vdot.
        the fist half of the list is xdot, the second half is vdot (for us vdot = g)
        x = [positions, velocities]
        xdot = [velocities, g]
        With time, this cuold be improved (atleast for clairity)
        adapted from https://mchouza.wordpress.com/2011/11/14/the-fall-of-the-slinky-i/
        '''
        xdot = [x[self.NX + i] if i < self.NX else self.G for i in range(2 * self.NX)]
        for i in range(self.NX - 1):
            a = self.sk * (x[i + 1] - x[i] - self.sl) / (self.mm * self.G)
            xdot[self.NX + i] += a
            xdot[self.NX + i + 1] -= a
        if self.held:
            xdot[self.NX] = 0.0 
        return xdot

    def rk4_step(self, x):
        # This operates on the whole list of positions at time t_i and returns the list of positions at time t_{i+1}
        k1 = self.get_xdot(x)
        k2 = self.get_xdot([x[i] + self.dt * k1[i] / 2.0 for i in range(len(x))])
        k3 = self.get_xdot([x[i] + self.dt * k2[i] / 2.0 for i in range(len(x))])
        k4 = self.get_xdot([x[i] + self.dt * k3[i] for i in range(len(x))])
        return [x[i] + self.dt * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0 for i in range(len(x))]
    
    def get_COM(self, x):
        return sum(x[:self.NX])/self.NX