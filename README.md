# SlugCode2

SlugCode2 is in early development stage.

## Install

```sh
mkdir build
cd build
../utils/setupSlug.py --fdm2d
# edit Makefile.header
make
```

It will create 2D/Hydro/FDM solver.
If you want to have 3D solver, use `--fdm3d` flag instead.

## Run

`SlugCode2` takes initial parameters from an external text file.
For example, if you'd like to use `vortex.init` file containing *all* initial parameters;

```sh
mpirun -np 4 ./slugCode2 ./vortex.init
```

The default filename for initial parameters is `slug.init`

```sh
mpirun -np 4 ./slugCode2   # this will read `slug.init`
```

`setupSlug.py` will copy every example `.init` files under `/inits` directory,
including `slug.init` under `src/{2,3}D`.

The variable name and its value in `.init` file should be separated with **single space**.

## Plot

I recommend to use [HDF5](https://www.hdfgroup.org) file output with parallel I/O.
For instance, you may specify initial parameters as:
```
...
# IO type
sim_hdf5 .true.
sim_pIO .true.
...
```
then, `SlugCode2` will write HDF5 format output in parallel.
This requires the HDF5 library with parallel support.

In order to read the data, you can use a dedicated Python module, [`slugger`](python/README.md).
This module reads the data from `SlugCode2` and returns `numpy` array:
```Python console
>>> import slugger as slug
>>> d = slug.load_data2d('vortex_10001.slug')
>>> d.dens
array([[0.9998759 , 0.99978231, 0.99969467, ..., 0.99979528, 0.99989642,
        0.99992699],
       [0.99992457, 0.99983486, 0.99973422, ..., 0.99981307, 0.99992151,
        0.99996232],
       [0.99992084, 0.99988414, 0.99981525, ..., 0.99978142, 0.99985646,
        0.9999095 ],
       ...,
       [0.99956576, 0.99961543, 0.9996903 , ..., 0.99964084, 0.99957778,
        0.99955136],
       [0.99966505, 0.99966503, 0.99968608, ..., 0.99968259, 0.99967797,
        0.99967071],
       [0.99978019, 0.99972438, 0.99968363, ..., 0.99974389, 0.99980183,
        0.99981395]])
```

`plotter/bin/plot{2,3}d` scripts use the `slugger` module to plot the data:
```sh
plotter/bin/plot2d cont vortex_10001.slug   # this will print density contour map
```

Also, you can use [visit](https://wci.llnl.gov/simulation/computer-codes/visit)
for more serious plottings.
Make sure you installed [SlugCode2-extension](plotter/visit/SlugCode/README.md) first.

## Supported schemes

### Finite Difference Method (FDM)
The conventional FDM scheme with Rusanov Lax-Friedrichs flux splitting and flux reconstruction.
See details on [A. Mignone *et. al.* Journal of Computational Physics 229.17 (2010): 5896-5920.][fdm]
 - Dimension: 2D and 3D
 - Spatial method: WENO-{JS, Z}, [GP-WENO][gp-weno]
 - Temporal method: Runge-Kutta {2, 3, 4}, [SF-PIF{3, 4}][sfpif]


[fdm]: https://doi.org/10.1016/j.jcp.2010.04.013
[gp-weno]: https://doi.org/10.1016/j.jcp.2018.12.028
[sfpif]: https://arxiv.org/abs/2006.00096
-------------------------------------------------------------------

## TODOs

- [ ] Implement internal boundary: `BDRY_VAR`
- [x] Add example `init` files.
- [x] Add abort conditions for `gr_ngc`
    - [x] `if gr_ngc > [gr_nx, gr_ny]: abort`
    - [x] `if gr_ngc < # of guad cell needed: abort`
- [ ] Rearrange folder structure -> eg:`FDM/[spatial temporal]`
- [ ] Makefile in each directory to handle `OBJS` accordingly
    - [ ] Edit `setup.py` to write relevant `Makefile`

### Do FVM

#### Spatial

- [ ] WENO
- [ ] GP-WENO

#### Temporal

- [ ] RK2
- [ ] RK3
- [ ] RK4
- [ ] PIF-SF
