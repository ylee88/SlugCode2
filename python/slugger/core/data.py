""" plotting tools for slugCode """
from copy import copy
from mpl_toolkits import axes_grid1  # for colorbar
import numpy as np
import matplotlib.pyplot as plt
import h5py as hdf


def get_data(data, var):
    """Return the data of given variable name: [dens, pres, velx, vely]"""
    try:
        var_data = getattr(data, var)
    except AttributeError:
        print("Error: " + var + " not found")
        raise

    return var_data


def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """draw color bar in proper position and size"""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1 / aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)

    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


class data1d:

    class ascii:

        def __init__(self, file_name):
            self.filename = file_name
            self.raw = np.loadtxt(file_name)
            self.x = self.raw[:, 0]
            self.dens = self.raw[:, 1]
            self.velx = self.raw[:, 2]
            self.pres = self.raw[:, 3]
            self.eint = self.raw[:, 4]

        def plot(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            ax.plot(self.x, var_data, label=self.filename, *args, **kwargs)


def load_data2d(file_name):

    if file_name.endswith('.dat'):
        return data2d.ascii(file_name)
    elif file_name.endswith('.slug'):
        return data2d.hdf5(file_name)
    else:
        raise Exception("Unrecognized data file extension.")


def load_data3d(file_name):

    if file_name.endswith('.dat'):
        return data3d.ascii(file_name)
    elif file_name.endswith('.slug'):
        return data3d.hdf5(file_name)
    else:
        raise Exception("Unrecognized data file extension.")


class data2d:

    class ascii:

        def __init__(self, file_name):
            self.raw = np.loadtxt(file_name)
            # self.bins = int(np.sqrt(np.size(self.raw[:, 0])))
            self.xbins = int(np.shape(np.unique(self.raw[:, 0]))[0])
            self.ybins = int(np.shape(np.unique(self.raw[:, 1]))[0])
            self.x = np.reshape(self.raw[:, 0], (self.xbins, self.ybins)).T
            self.y = np.reshape(self.raw[:, 1], (self.xbins, self.ybins)).T
            self.dens = np.reshape(self.raw[:, 2], (self.xbins, self.ybins)).T
            self.velx = np.reshape(self.raw[:, 3], (self.xbins, self.ybins)).T
            self.vely = np.reshape(self.raw[:, 4], (self.xbins, self.ybins)).T
            self.pres = np.reshape(self.raw[:, 5], (self.xbins, self.ybins)).T

        def __add__(self, other):
            total = copy(self)
            total.raw = 0.
            # total.bins = self.bins
            total.xbins = self.xbins
            total.ybins = self.ybins
            total.x = self.x
            total.y = self.y

            total.dens = self.dens + other.dens
            total.velx = self.velx + other.velx
            total.vely = self.vely + other.vely
            total.pres = self.pres + other.pres

            return(total)

        def __sub__(self, other):
            total = copy(self)
            total.raw = 0.
            # total.bins = self.bins
            total.xbins = self.xbins
            total.ybins = self.ybins
            total.x = self.x
            total.y = self.y

            total.dens = self.dens - other.dens
            total.velx = self.velx - other.velx
            total.vely = self.vely - other.vely
            total.pres = self.pres - other.pres

            return(total)

        def edge_grid(self):

            dx = self.x[0, 1] - self.x[0, 0]
            dy = self.y[1, 0] - self.y[0, 0]

            # shift center points by delta/2
            # this is one dimensional array
            xx = self.x[0, :] - dx/2.
            yy = self.y[:, 0] - dy/2.

            # append last point
            xi = np.append(xx, xx[-1]+dx)
            yi = np.append(yy, yy[-1]+dy)

            return xi, yi

        def check_symmetric(self, var, tol=1e-8):

            var_data = get_data(self, var)

            if np.shape(var_data) == np.shape(var_data.T):
                sym_diag = np.allclose(var_data, var_data.T, atol=tol)
                sym_lr = np.allclose(var_data, np.fliplr(var_data), atol=tol)
                sym_ud = np.allclose(var_data, np.flipud(var_data), atol=tol)
                symmetric = sym_lr and sym_ud and sym_diag
                print(sym_lr, sym_ud, sym_diag, symmetric)

                if(not symmetric):
                    print("[[[ WARN: DATA IS ASYMMETRIC ]]]")

                return(symmetric)
            else:
                return ""

        def plot_cmap(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            xi, yi = self.edge_grid()

            im = ax.pcolormesh(xi, yi, var_data, *args, **kwargs)
            ax.set_title(var)
            add_colorbar(im)

        def plot_contour(self, var, ax, nlevel=20, *args, **kwargs):

            var_data = get_data(self, var)

            extent = (np.amin(self.x), np.amax(self.x),
                      np.amin(self.y), np.amax(self.y))

            levels = np.linspace(var_data.min(), var_data.max(), nlevel)
            # levels = np.log(levels)
            # levels = np.logspace(var_data.min(), var_data.max(), nlevel)

            ax.contour(var_data, extent=extent, levels=levels, *args, **kwargs)
            ax.set_title(var)

        def plot_slicex(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.ybins % 2 == 0:
                pos1 = int(self.ybins/2)
                pos2 = int(self.ybins/2) + 1
                slice_data = 0.5*(var_data[pos1, :] + var_data[pos2, :])
            else:
                pos = int(self.ybins/2) + 1
                slice_data = var_data[pos, :]

            ax.plot(self.x[0, :], slice_data, *args, **kwargs)
            ax.set_title(var + ":slice in x")

        def plot_slicey(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.xbins % 2 == 0:
                pos1 = int(self.xbins/2)
                pos2 = int(self.xbins/2) + 1
                slice_data = 0.5*(var_data[:, pos1] + var_data[:, pos2])
            else:
                pos = int(self.xbins/2) + 1
                slice_data = var_data[:, pos]

            ax.plot(self.y[:, 0], slice_data, *args, **kwargs)
            ax.set_title(var + ":slice in y")

        def plot_slice45(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.xbins == self.ybins:
                slice45 = np.zeros(self.xbins)
                for i in range(self.xbins):
                    slice45[i] = var_data[i, i]

                ax.plot(np.sqrt(2)*self.x, slice45, *args, **kwargs)
                ax.set_title(var + ":slice in 45 degree")
            else:
                return

        def plot_symlr(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)
            data_lr = np.fliplr(var_data)

            diff = np.abs(var_data - data_lr)
            xi, yi = self.edge_grid()

            im = ax.pcolormesh(xi, yi, diff, *args, **kwargs)
            ax.set_title('diff_LR : ' + var + '\nsum : ' +
                         format(np.sum(diff), 'e'))
            add_colorbar(im)

        def plot_symud(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)
            data_ud = np.flipud(var_data)

            diff = np.abs(var_data - data_ud)
            xi, yi = self.edge_grid()

            im = ax.pcolormesh(xi, yi, diff, *args, **kwargs)
            ax.set_title('diff_UD : ' + var + '\nsum : ' +
                         format(np.sum(diff), 'e'))
            add_colorbar(im)

        def plot_sym45(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)
            data_45 = var_data.T

            diff = np.abs(var_data - data_45)
            xi, yi = self.edge_grid()

            im = ax.pcolormesh(xi, yi, diff, *args, **kwargs)
            ax.set_title('diff_45 : ' + var + '\nsum : ' +
                         format(np.sum(diff), 'e'))
            add_colorbar(im)

    class hdf5(ascii):

        def __init__(self, file_name):

            f = hdf.File(file_name, "r")
            self.raw = np.array(f.get('prim_vars'))
            self.ybins, self.xbins = np.shape(self.raw[:, :, 0])
            self.x = f.get('xCoord')
            self.y = f.get('yCoord')

            self.dens = self.raw[:, :, 0]
            self.velx = self.raw[:, :, 1]
            self.vely = self.raw[:, :, 2]
            self.pres = self.raw[:, :, 3]

            self.eTime = np.array(f.get('eTime'))[0]

# below functions should be redefined as the dimension of xCoord and yCoord are
# different with ascii data's -> TODO: fixit!
        def edge_grid(self):

            dx = self.x[1] - self.x[0]
            dy = self.y[1] - self.y[0]

            # shift center points by delta/2
            # this is one dimensional array
            xx = self.x[:] - dx/2.
            yy = self.y[:] - dy/2.

            # append last point
            xi = np.append(xx, xx[-1]+dx)
            yi = np.append(yy, yy[-1]+dy)

            return xi, yi

        def plot_slicex(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.ybins % 2 == 0:
                pos1 = int(self.ybins/2)
                pos2 = int(self.ybins/2) + 1
                slice_data = 0.5*(var_data[pos1, :] + var_data[pos2, :])
            else:
                pos = int(self.ybins/2) + 1
                slice_data = var_data[pos, :]

            ax.plot(self.x, slice_data, *args, **kwargs)
            ax.set_title(var + ":slice in x")

        def plot_slicey(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.xbins % 2 == 0:
                pos1 = int(self.xbins/2)
                pos2 = int(self.xbins/2) + 1
                slice_data = 0.5*(var_data[:, pos1] + var_data[:, pos2])
            else:
                pos = int(self.xbins/2) + 1
                slice_data = var_data[:, pos]

            ax.plot(self.y, slice_data, *args, **kwargs)
            ax.set_title(var + ":slice in y")

        def plot_slice45(self, var, ax, *args, **kwargs):

            var_data = get_data(self, var)

            if self.xbins == self.ybins:
                slice45 = np.zeros(self.xbins)
                for i in range(self.xbins):
                    slice45[i] = var_data[i, i]

                ax.plot(np.sqrt(2)*self.x, slice45, *args, **kwargs)
                ax.set_title(var + ":slice in 45 degree")
            else:
                return


class data3d:

    class ascii:

        def __init__(self, file_name):
            self.raw = np.loadtxt(file_name)

            self.xbins = int(np.shape(np.unique(self.raw[:, 0]))[0])
            self.ybins = int(np.shape(np.unique(self.raw[:, 1]))[0])
            self.zbins = int(np.shape(np.unique(self.raw[:, 2]))[0])
            self.x = np.reshape(self.raw[:, 0], (self.xbins, self.ybins, self.zbins)).T
            self.y = np.reshape(self.raw[:, 1], (self.xbins, self.ybins, self.zbins)).T
            self.z = np.reshape(self.raw[:, 2], (self.xbins, self.ybins, self.zbins)).T
            self.dens = np.reshape(self.raw[:, 3], (self.xbins, self.ybins, self.zbins)).T
            self.velx = np.reshape(self.raw[:, 4], (self.xbins, self.ybins, self.zbins)).T
            self.vely = np.reshape(self.raw[:, 5], (self.xbins, self.ybins, self.zbins)).T
            self.velz = np.reshape(self.raw[:, 6], (self.xbins, self.ybins, self.zbins)).T
            self.pres = np.reshape(self.raw[:, 7], (self.xbins, self.ybins, self.zbins)).T

        def __add__(self, other):
            total = copy(self)
            total.raw = 0.
            total.xbins = self.xbins
            total.ybins = self.ybins
            total.zbins = self.zbins
            total.x = self.x
            total.y = self.y
            total.z = self.z

            total.dens = self.dens + other.dens
            total.velx = self.velx + other.velx
            total.vely = self.vely + other.vely
            total.velz = self.velz + other.velz
            total.pres = self.pres + other.pres

            return(total)

        def __sub__(self, other):
            total = copy(self)
            total.raw = 0.
            total.xbins = self.xbins
            total.ybins = self.ybins
            total.ybins = self.zbins
            total.x = self.x
            total.y = self.y
            total.z = self.z

            total.dens = self.dens - other.dens
            total.velx = self.velx - other.velx
            total.vely = self.vely - other.vely
            total.velz = self.velz - other.velz
            total.pres = self.pres - other.pres

            return(total)

    class hdf5(ascii):

        def __init__(self, file_name):

            f = hdf.File(file_name, "r")
            self.raw = np.array(f.get('prim_vars'))
            self.zbins, self.ybins, self.xbins = np.shape(self.raw[:, :, :, 0])
            self.x = f.get('xCoord')
            self.y = f.get('yCoord')
            self.z = f.get('zCoord')

            self.dens = self.raw[:, :, :, 0].T
            self.velx = self.raw[:, :, :, 1].T
            self.vely = self.raw[:, :, :, 2].T
            self.velz = self.raw[:, :, :, 3].T
            self.pres = self.raw[:, :, :, 4].T

            self.eTime = np.array(f.get('eTime'))[0]
