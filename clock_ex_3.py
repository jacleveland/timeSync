import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from bisect import bisect_left

"""
    Simplified LCRD
    10 --> 1:[2, 3, 4] --> 5 --> 100
"""

class clock:
    def __init__(self, time_start, dt):
        self.time = time_start
        self.dt = dt
    
    def increment_clock(self):
        self.time += self.dt
    
    def set_time(self, time):
        self.time = time
    
    def set_dt(self, dt):
        self.dt = dt

class node:
    def __init__(self, id, time_data, send_sched, update_sched, coef):
        self.id = id
        self.clock = clock(time_data[0], time_data[1])
        self.send_sched = send_sched
        self.update_sched = update_sched
        self.coef = coef
        self.data = {}
        self.nbhd = []

    # All neighborhood nodes (need to be added after initialization)
    def add_nbhd(self, nbhd):
        self.nbhd = nbhd
    
    # To string
    def __repr__(self):
        st = f"Node id: {str(self.id)} Current Time: {self.clock.time:.4f}"
        return st 

    # Send clock data
    def send(self, dest):
        dest.recieve(self.id, self.clock.time)

    # Recieve clock data
    def recieve(self, other_id, other_time):
        self.data[other_id] = other_time

    # Check if a time is in the schedule (gives more problems than you'd hope)
    def validTime(self, time, sched):
        idx = bisect_left(sched, time)
        if idx < len(sched) and np.abs(sched[idx] - time) < 1e-4:
            return True
        
        # This might be off...
        if (time >= 40023.050 and time <= 40023.090):
            print("BUG WITH VERIFIER", self.id, time)

        return False

    def schedule_send(self, conn):
        if self.validTime(self.clock.time, self.send_sched):
            for nbhs in self.nbhd:
                if nbhs in conn:
                    self.send(nbhs)


    def update(self):
        if self.validTime(self.clock.time, self.send_sched):
            self.euler()
        
        self.clock.increment_clock()

    def euler(self):
        upd = self.coef * len(self.data) * self.clock.time
        for i in self.data:
            upd -= self.coef * self.data[i]
        
        self.clock.set_time(self.clock.time - upd)
        self.data = {}


class lcrd:
    def __init__(self, ext_nodes, sched):
        self.id = 1
        # Create nodes
        # node id, [start time, increment value], send schedule, update schedule, coefficient
        n2 = node(2, [40023.043, 0.0001], sched, sched, 0.487)
        n3 = node(3, [40023.057, 0.0001], sched, sched, 0.487)
        n4 = node(4, [40023.064, 0.0001], sched, sched, 0.487)
        n5 = node(5, [40023.044, 0.0001], sched, sched, 0.31)
        
        # Add connections within nodes
        n2.add_nbhd([ext_nodes[0], n5])
        n3.add_nbhd([ext_nodes[0], n5])
        n4.add_nbhd([ext_nodes[0], n5])
        n5.add_nbhd([n2, n3, n4, ext_nodes[1]])

        self.nodes = [n2, n3, n4, n5]
        
        # Currently connected to ground station 
        self.gs = node(-1, [-1, -1], -1, -1, -1)
        
    def recieve(self, other_id, other_time):
        if (other_id == 10):
            self.gs.recieve(other_id, other_time)
        elif(other_id == 100):
            # node 5
            self.nodes[3].recieve(other_id, other_time)

    # Represent switching between operational satellite by using computer time
    def updateStation(self, curr):
        if (curr < 20 and self.gs.id != 2):
            self.gs = self.nodes[0]
            print("Current Ground Station Available:", 2)
        elif (curr >= 20 and curr < 40 and self.gs.id != 3):
            self.gs = self.nodes[1]
            print("Current Ground Station Available:", 3)
        elif (curr >= 40 and self.gs.id != 4):
            self.gs = self.nodes[2]
            print("Current Ground Station Available:", 4)    
        
        self.gs.coef = .19


    def run(self, time):
        # Switch ground station if needed
        self.updateStation(time)

        # Transmit data within and outside lcrd
        for node in self.nodes:
            if node.id != self.gs.id and node.id != 5:
                node.schedule_send([node.nbhd[1]])
            else:
                node.schedule_send(node.nbhd)
        


class network:
    def __init__(self):
        # node id, [start time, increment value], send schedule, update schedule, coefficient
        sched = np.arange(40023.050, 40023.0901, 0.0001)
        self.n10 = node(10, [40023.074, 0.0001], sched, sched, 0.487)
        self.n100 = node(100, [40023.055, 0.0001], sched, sched, 0.487)
        self.lcrd = lcrd([self.n10, self.n100], sched)
        
        self.n10.add_nbhd([self.lcrd])
        self.n100.add_nbhd([self.lcrd])
    
    # Set up lists for plotting results
    def plotter_setup(self, nodes):
        plots = []
        for i in range(len(nodes)):
            plots.append([nodes[i].clock.time])

        return plots
    
    # Update lists for plotting results
    def plotter_upd(self, plots):
        plots[0].append(self.n100.clock.time)
        plots[1].append(self.n10.clock.time)
        plots[2].append(self.lcrd.nodes[0].clock.time)
        plots[3].append(self.lcrd.nodes[1].clock.time)
        plots[4].append(self.lcrd.nodes[2].clock.time)
        plots[5].append(self.lcrd.nodes[3].clock.time)

    # Plot time differences at each iteration
    def plotter_res(self, plots, c):
        # Plotter
        for i in range(len(c)):
            avg = np.average([plots[0][i], plots[1][i], plots[2][i], plots[3][i], plots[4][i], plots[5][i]])
            for j in range(len(plots)):
                plots[j][i] = plots[j][i] - avg

        plt.rcParams['font.size'] = 24
        plt.plot(c, plots[0])
        plt.plot(c, plots[1], color='maroon')
        plt.plot(c, plots[2], color='mediumpurple')
        plt.plot(c, plots[3], color='hotpink')
        plt.plot(c, plots[4], color='darkorange')
        plt.plot(c, plots[5], color='darkgreen')
        plt.legend(['n100', 'n10', 'n2', 'n3', 'n4', 'n5'])
        plt.xlabel("Iteration Count")
        plt.ylabel("Time Difference (Seconds)")
        plt.title("Difference from Mean Clock Time vs. Iteration Count (Ground Station 4)")
        plt.show()

    # Runs the simulation
    def run(self, dur):
        # Plotter
        c = [0]
        i = 0
        plots = self.plotter_setup([self.n100, self.n10, self.lcrd.nodes[0], self.lcrd.nodes[1], self.lcrd.nodes[2], self.lcrd.nodes[3]])
        start = timer()
        time = timer()
        while (np.abs(time - start) <= dur):
            i += 1
            c.append(i)
            lc_time = int(time) % 60
            self.lcrd.run(lc_time)
            
            # Send data
            self.n10.schedule_send([self.lcrd])
            self.n100.schedule_send([self.lcrd])

            # Update clocks
            self.n10.update()
            self.n100.update()
            for node in self.lcrd.nodes:
                node.update()

            self.plotter_upd(plots)
            time = timer()

        self.plotter_res(plots, c)
            
    # To String for network        
    def __repr__(self):
        string = "Network Data:\n" + repr(self.n10)
        for i in range(len(self.lcrd.nodes)):
            string += "\n" + repr(self.lcrd.nodes[i])
        string += "\n" + repr(self.n100)
        return string
        

n = network()
print(n)
# Run for approximately .005 seconds
n.run(.005)
print(n)