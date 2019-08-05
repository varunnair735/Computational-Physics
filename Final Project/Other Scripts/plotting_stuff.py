import matplotlib.pyplot as plt

fig9, ax9 = plt.subplots(1, 1, figsize = (10, 5),dpi=240)

ax9.imshow(p3,origin='lower',cmap='Blues')
ax9.set_title("Pressure ($Pa$)")
ax9.set_xlabel("$x\ (m)$")
ax9.set_ylabel("$y\ (m)$")
plt.savefig("{:2.0f}_wing_pressure.png".format(N))
