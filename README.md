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

## TODOs

- [ ] Implement internal boundary: `BDRY_VAR`
- [x] Add example `init` files.
- [x] Add abort conditions for `gr_ngc`
    - [x] `if gr_ngc > [gr_nx, gr_ny]: abort`
    - [x] `if gr_ngc < # of guad cell needed: abort`
- [ ] Rearrange folder structure -> eg:`FDM/[spatial temporal]`
- [ ] Makefile in each directory to handle `OBJS` accordingly
    - [ ] Edit `setup.py` to write relevant `Makefile`

### FDM

- Spatial: WENO, gp-WENO
- Temporal: RK[2, 3, 4], sfPIF[3, 4]

- [x] 2D
- [x] 3D

### FVM

#### Spatial

- [ ] WENO
- [ ] GP-WENO

#### Temporal

- [ ] RK2
- [ ] RK3
- [ ] RK4
- [ ] PIF-SF
