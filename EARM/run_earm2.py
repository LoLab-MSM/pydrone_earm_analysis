import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
from earm2_flat import model

tspan = np.linspace(0, 20000, 100)

sim = ScipyOdeSimulator(model, tspan=tspan).run().all

plt.plot(tspan, sim['cPARP'], label='cPARP')
plt.xlabel('Time(s)')
plt.ylabel('Population')
plt.legend()
plt.show()