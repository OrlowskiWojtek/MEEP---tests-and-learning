import meep as mp
import math
import numpy as np
import matplotlib.pyplot as plt

def main(n):
    # Some parameters to describe the geometry:
    eps = 12  # dielectric constant of waveguide
    kx = 0
    # The cell dimensions
    dpml = 20   # PML thickness (y direction only!)
    md = n*(2+1/math.sqrt(3))
    ycell=md+10 
    sy = 2*dpml + ycell # size of cell in y direction (perpendicular to wvg.)

    cell = mp.Vector3(0, sy)

    fcen = 0.125  # pulse center frequency
    df = 0.2

    wvlngth = 8
    wvlwidth = 8

    s = mp.Source(
        src=mp.GaussianSource(frequency=fcen, fwidth=df),
        component=mp.Ez,
        center=mp.Vector3(0,-ycell/2)
    )

    sim = mp.Simulation(
        cell_size=cell,
        geometry=None,
        sources=[s],
        k_point = mp.Vector3(kx),
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        resolution=30
    )

    nfreq = 1000  # number of frequencies at which to compute flux

    refl_fr = mp.FluxRegion(center=mp.Vector3(0,-ycell/2+1,0), direction=mp.Y)
    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    tran_fr = mp.FluxRegion(center=mp.Vector3(0,ycell/2-1,0), direction=mp.Y)
    tran = sim.add_flux(fcen, df, nfreq, tran_fr)

    pt = mp.Vector3(0,ycell/2-1)

    #sim.run(mp.to_appended("ezSim1", mp.at_every(0.6, mp.output_efield_z)), until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-3))

    # for normalization run, save flux fields data for reflection plane
    straight_refl_data = sim.get_flux_data(refl)

    straight_tran_flux = mp.get_fluxes(tran)

    sim.reset_meep()

    geometry = []

    for i in range(n):
        geometry.append(mp.Block(size = mp.Vector3(mp.inf,1/(math.sqrt(3)),mp.inf),center= mp.Vector3(0,-md/2+1/(2*math.sqrt(3))+i*(1/math.sqrt(3) + 2),0),
                                 material = mp.Medium(epsilon=eps)))
    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        sources=[s],
        k_point = mp.Vector3(kx),
        boundary_layers=[mp.PML(dpml, direction=mp.Y)],
        resolution=30
    )

    refl = sim.add_flux(fcen, df, nfreq, refl_fr)

    tran_fr = mp.FluxRegion(center=mp.Vector3(0,ycell/2-1,0), direction=mp.Y)
    tran = sim.add_flux(fcen, df, nfreq, tran_fr)

    sim.load_minus_flux_data(refl, straight_refl_data)

    #sim.run(mp.at_beginning(mp.output_epsilon),mp.to_appended("ezSim2", mp.at_every(0.6, mp.output_efield_z)),until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, pt, 1e-3))

    bend_refl_flux = mp.get_fluxes(refl)
    bend_tran_flux = mp.get_fluxes(tran)

    flux_freqs = mp.get_flux_freqs(refl)

    fr = []
    Rs = []
    Ts = []
    for i in range(nfreq):
        fr = np.append(fr, flux_freqs[i])
        Rs = np.append(Rs,-bend_refl_flux[i]/straight_tran_flux[i])
        Ts = np.append(Ts,bend_tran_flux[i]/straight_tran_flux[i])

    if mp.am_master():
        plt.figure()
        plt.title("n="+str(n))
        #plt.plot(fr,Rs,'bo-',label='reflectance')
        plt.plot(fr,Ts,'ro-',label='transmittance')
        #plt.plot(fr,1-Rs-Ts,'go-',label='loss')
        plt.axis([fcen - df/2-0.05 , fcen + df/2+0.05, 0, 1])
        plt.xlabel("Frequency (1/μm)")
        plt.legend(loc="upper right")

        plt.savefig("gallery/vnTRANS"+str(n)+".png")

if __name__ == "__main__":
    for i in range (0,20):
        main(i)