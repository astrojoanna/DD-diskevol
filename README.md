This code was written by C.P. Dullemond and J. Drążkowska.

The disk evolution part was used in:

* [Dullemond, Natta & Testi (2006) ApJ 645, L69](https://ui.adsabs.harvard.edu/abs/2006ApJ...645L..69D/abstract)
* [Dullemond, Apai & Walch (2006) ApJ 640, L67](https://ui.adsabs.harvard.edu/abs/2006ApJ...640L..67D/abstract)

The disk evolution + dust evolution modules were used in:

* [Drążkowska & Dullemond (2018) A&A 614, A62](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..62D/abstract)
* [Lichtenberg et al. (2021) Science 371, 6527, 365](https://science.sciencemag.org/content/371/6527/365)

HOW TO RUN THE CODE

Review `Makefile` to include your Fortran compiler and compiler flags.

Compile with `make`.

Run with `./diskevol params.par`, where `params.par` is parameter file.
The parameter file included in this repository corresponds to fiducial setup from [Drążkowska & Dullemond (2018) A&A 614, A62](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..62D/abstract), which was used in Lichtenberg et al. (2020).

OUTPUT

All outputs are in CGS units.
One file includes data corresponding to one variable:

* time.info         number of outputs produced,
* grid.info			    radial distance to the central star,
* sigma.dat			    gas surface density,
* sigmad.dat		    dust surface density,
* sigmaice.dat		  ice component surface density,
* sigmavap.dat		  water vapour surface density,
* sigmaplts.dat		  planetesimals surface density,
* dustsize.dat		  representative dust grain size,
* stokesnr.dat		  representative Stokes number,
* stfrag.dat			  turbulent fragmentation limited Stokes number,
* stdf.dat 			    drift fragmentation limited Stokes number,
* stdrift.dat 		  drift barrier limited Stokes number,
* etamid.dat			  mid plane dust-to-gas ratio,
* vdust.dat			    dust velocity,
* mflux.dat			    dust mass flux through cell interfaces,
* vdrift.dat        maximum velocity of radial drift = headwind speed,
* velo.dat          radial velocity of gas,
* temperature.dat 	midplane temperature,
* visc.dat			    gas viscosity,
* alpha.dat			    gas turbulence strength parameter,
* mdot.dat			    gas mass flow through cell interface,
* mstar.dat         mass of the central star,
* time.dat		      output time,

For time dependent variables, the files are appended each time the output is produced, the outputs are separated by blank line. One row corresponds to one radial cell.

LICENSE

This code is licensed under the GNU General Public License v3. Among other things, this means that if you modified this code for a publication, the original code needs to be cited (see below), the changes need to be stated, and the new code needs to be open sourced under the same license.

CREDITS

If you make use of this code, or any of its parts, please cite this repository. 
If you use the disk evolution part, please cite [Dullemond, Natta & Testi (2006) ApJ 645, L69](https://ui.adsabs.harvard.edu/abs/2006ApJ...645L..69D/abstract).
If you use the dust evolution part, please cite [Drążkowska & Dullemond (2018) A&A 614, A62](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..62D/abstract). 
