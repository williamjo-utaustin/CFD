import numpy as np
import matplotlib.pyplot as plt


y = np.linspace(0,0.1,100)
u = 1.5*(1 - (y-1)**2/(0.1)**2)

plt.plot(y,u)
plt.show()



