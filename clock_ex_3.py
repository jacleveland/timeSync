import simpy
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import time as tm
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
        return "Node id: " + str(self.id) + " Current Time: " + str(self.clock.time)

    # Send clock data
    def send(self, dest):
        #TODO: CHECK IF SCHEDULE SAYS TO SEND DATA TO OTHER NEIGHBORS
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
        #print(np.any(np.abs(np.full_like(self.send_sched, self.clock.time) - self.send_sched) <= 1e-6))
        if self.validTime(self.clock.time, self.send_sched):
            #print("Send:", self.id, "at time", self.clock.time, end=" ")
            for nbhs in self.nbhd:
                if nbhs in conn:
                    #print("to:", nbhs.id, end=" ")
                    self.send(nbhs)
            #print()
        #else:
            #print(self.clock.time)

    def update(self):
        if self.validTime(self.clock.time, self.send_sched):
            self.euler()
        
        self.clock.increment_clock()

    def euler(self):
        #print(self.id, len(self.data), self.data)
        upd = self.coef * len(self.data) * self.clock.time
        for i in self.data:
            upd -= self.coef * self.data[i]
        
        self.clock.set_time(self.clock.time - upd)
        self.data = {}


class lcrd:
    def __init__(self, ext_nodes, sched):
        self.id = 1
        # Create nodes
        #sched = np.arange(40023.044, 40024.120, .001)
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
        #print("lcrd recieve", other_id)
        if (other_id == 10):
            #print(self.gs.id)
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
    

    def run(self):
        start = timer()
        time = timer()
        #for i in range(38):
        # Run for approximately .1 second
        while (time - start < .1): 
            #print(time)
            #st = timer()
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
            #print(timer() - st)
            #self.n100.euler()
            #self.n100.clock.increment_clock()
            #self.n10.euler()
            #self.n10.clock.increment_clock()
            #for node in self.lcrd.nodes:
                #node.euler()
                #node.clock.increment_clock()
            tm.sleep(.0001)
            time = int(timer())
            
    # To String for network        
    def __repr__(self):
        string = "Network Data:\n" + repr(self.n10)
        for i in range(len(self.lcrd.nodes)):
            string += "\n" + repr(self.lcrd.nodes[i])
        string += "\n" + repr(self.n100)
        return string
        


## NO STRATUM
"""
n10 = node(10, [40023.054, 0.001], -1, -1, 0.5 * np.array([1, -1]))
n2 = node(2, [40023.043, 0.001], -1, -1, 0.4999 * np.array([[-1, 2, -1], [0, 1, -1]]))
n3 = node(3, [40023.057, 0.001], -1, -1, 0.4999 * np.array([[-1, 2, -1], [0, 1, -1]]))
n4 = node(4, [40023.064, 0.001], -1, -1, 0.4999 * np.array([[-1, 2, -1], [0, 1, -1]]))
n5 = node(5, [40023.044, 0.001], -1, -1, 0.31 * np.array([-1, -1, -1, 4, -1]))
n100 = node(100, [40023.055, 0.001], -1, -1, 0.5 * np.array([-1, 1]))

print(n10)
"""

n = network()
print(n)
n.run()
print(n)


"""
# Represent switching between operational satellite by using computer time
def updateStation(time):
    gs = -1
    curr = int(time) % 60
    #curr = int(timer()) % 60
    print(curr)
    if (curr < 20 and gs != 2):
        gs = 2
        print("Current Ground Station Available:", gs)
    elif (curr >= 20 and curr < 40 and gs != 3):
        gs = 3
        print("Current Ground Station Available:", gs)
    elif (curr >= 40 and gs != 4):
        gs = 4
        print("Current Ground Station Available:", gs)    
    #time.sleep(1)
"""


"""
# Use this iteration to increment clocks
# change dt for different rates
for i in range(1):
    #start = timer()
    # around 2e-5 to 3e-5 to increment all the clocks
    # So make sure you sleep for a value greater than that
    n10.clock.increment_clock()
    n2.clock.increment_clock()
    n3.clock.increment_clock()
    n4.clock.increment_clock()
    n5.clock.increment_clock()
    n100.clock.increment_clock()
    #end = timer()
    #print(end - start)
    # .001 sec / sec
    time.sleep(1)
"""
