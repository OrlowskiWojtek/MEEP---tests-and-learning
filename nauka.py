import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import os

os.system("mkdir gifs")

resolution = 16 # pixels/um

sx = 32
sy = 32

dpml = 4.0

cell = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(dpml)]

geometry = []

symmetries = [mp.Mirror(mp.Y)]

fcen = 1/8 # pulse center frequency
df = 0.1   # pulse width (in frequency)
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df,is_integrated=True),
                     component=mp.Ez,
                     center=mp.Vector3(-sx/2 + dpml + 1,0,0),
                     size=mp.Vector3(0,sy-2,0))]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=None,
                    sources=sources,
                    resolution=resolution)


nfreq = 100  # number of frequencies at which to compute flux

refl_fr = mp.FluxRegion(center=mp.Vector3(-sx/2 + dpml + 2,0,0), size=mp.Vector3(0,sy-2*dpml,0))
refl = sim.add_flux(fcen, df, nfreq, refl_fr)

tran_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - dpml - 2,0,0), size=mp.Vector3(0,sy-2*dpml,0))
tran = sim.add_flux(fcen, df, nfreq, tran_fr)

pt = mp.Vector3(sx/2-dpml-1,0)

sim.run(mp.to_appended("ezSim1", mp.at_every(0.6, mp.output_efield_z)),until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

# for normalization run, save flux fields data for reflection plane
straight_refl_data = sim.get_flux_data(refl)

straight_tran_flux = mp.get_fluxes(tran)

sim.reset_meep()

geometry = []

for i in range(3):
    geometry.append(mp.Block(mp.Vector3(2/(np.sqrt(12)),mp.inf,mp.inf), center = mp.Vector3(2*i,0,0), material = mp.Medium(epsilon=12)))


sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

# reflected flux
refl = sim.add_flux(fcen, df, nfreq, refl_fr)

tran_fr = mp.FluxRegion(center=mp.Vector3(sx/2-dpml-2,0,0), size=mp.Vector3(0,sy-2*dpml,0))
tran = sim.add_flux(fcen, df, nfreq, tran_fr)

# for normal run, load negated fields to subtract incident from refl. fields
sim.load_minus_flux_data(refl, straight_refl_data)

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ezSim2", mp.at_every(0.6, mp.output_efield_z)),until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

bend_refl_flux = mp.get_fluxes(refl)
bend_tran_flux = mp.get_fluxes(tran)

flux_freqs = mp.get_flux_freqs(refl)

wl = []
Rs = []
Ts = []
for i in range(nfreq):
    wl = np.append(wl, 1/flux_freqs[i])
    Rs = np.append(Rs,-bend_refl_flux[i]/straight_tran_flux[i])
    Ts = np.append(Ts,bend_tran_flux[i]/straight_tran_flux[i])

if mp.am_master():
    plt.figure()
    plt.plot(wl,Rs,'bo-',label='reflectance')
    plt.plot(wl,Ts,'ro-',label='transmittance')
    plt.plot(wl,1-Rs-Ts,'go-',label='loss')
    plt.axis([4, 20, 0, 1])
    plt.xlabel("wavelength (μm)")
    plt.legend(loc="upper right")
    plt.savefig('gallery/plot.png')