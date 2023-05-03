import meep as mp
import math
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Some parameters to describe the geometry:
    eps = 12  # dielectric constant of waveguide
    kx = 0
    # The cell dimensions
    sy = 20 # size of cell in y direction (perpendicular to wvg.)
    dpml = 1  # PML thickness (y direction only!)

    cell = mp.Vector3(0, sy)

    fcen = 0.1  # pulse center frequency
    df = 0.09  # pulse freq. width: large df = short impulse

    s = mp.Source(
        src=mp.GaussianSource(fcen, fwidth=df),
        component=mp.Ex,
        center=mp.Vector3(0,-7)
    )

    sim = mp.Simulation(
        cell_size=cell,
        geometry=None,
        sources=[s],
        k_point = mp.Vector3(kx),
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        resolution=20
    )

    nfreq = 100  # number of frequencies at which to compute flux

    refl_fr = mp.FluxRegion(center=mp.Vector3(0,-5,0), direction=mp.Y)
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    tran_fr = mp.FluxRegion(center=mp.Vector3(0,sy/2 - dpml - 2,0), direction=mp.Y)
    tran = sim.add_flux(fcen, df, nfreq, tran_fr)

    pt = mp.Vector3(0,sy/2-dpml-1)

    sim.run(mp.to_appended("ezSim1", mp.at_every(0.6, mp.output_efield_x)), until_after_sources=mp.stop_when_fields_decayed(50,mp.Ex,pt,1e-3))

    # for normalization run, save flux fields data for reflection plane
    straight_refl_data = sim.get_flux_data(refl)

    straight_tran_flux = mp.get_fluxes(tran)

    sim.reset_meep()

    geometry = []

    for i in range(3):
        geometry.append(mp.Block(size = mp.Vector3(mp.inf,2/(math.sqrt(12)),mp.inf),center= mp.Vector3(0,-3+2*i,0),
                                 material = mp.Medium(epsilon=eps)))
    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        sources=[s],
        k_point = mp.Vector3(kx),
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        resolution=20
    )

    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    tran_fr = mp.FluxRegion(center=mp.Vector3(0,sy/2 - dpml - 2,0), direction=mp.Y)
    tran = sim.add_flux(fcen, df, nfreq, tran_fr)

    sim.load_minus_flux_data(refl, straight_refl_data)

    sim.run(mp.at_beginning(mp.output_epsilon),mp.to_appended("ezSim2", mp.at_every(0.6, mp.output_efield_x)),until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, pt, 1e-3))

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
        plt.show()
    
    omega = 0.125

    s = [
            mp.Source(
            src=mp.GaussianSource(omega, fwidth=0.01),
            component=mp.Hz,
            center=mp.Vector3(0, -5),
        )
    ]
    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        sources=s,
        k_point=mp.Vector3(kx),
        resolution=0,
    )
    f = plt.figure(dpi=100)
    animate = mp.Animate2D(fields=mp.Hx, f=f, normalize=True, realtime=False)
    sim.run(mp.at_every(5, animate), until_after_sources=1)
    animate.to_mp4(10, "test.mp4")
    plt.close()

if __name__ == "__main__":
    main()