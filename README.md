# SlugCode2

SlugCode2 is in early development stage.

## Install

```sh
mkdir build
cd build
../utils/setupSlug.py
# edit Makefile.header
make
```

It will create 2D/Hydro/FDM solver.

## TODOs

- [ ] Implement internal boundary: `BDRY_VAR`
- [ ] Add example `init` files.
- [ ] Add abort conditions for `gr_ngc`
    - [ ] `if gr_ngc > [gr_nx, gr_ny]: abort`
    - [ ] `if gr_ngc < # of guad cell needed: abort`
- [ ] Rearrange folder structure -> eg:`FDM/[spatial temporal]`
- [ ] Makefile in each directory to handle `OBJS` accordingly
    - [ ] Edit `setup.py` to write relevant `Makefile`

### FDM

#### Spatial

- [x] WENO
- [ ] GP

#### Temporal

- [x] RK2
- [x] RK3
- [x] RK4
- [x] PIF-SF

### FVM
