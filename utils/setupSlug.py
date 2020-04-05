#!/usr/bin/env python3
import os
import shutil


class slug_config:

    class item:

        def __init__(self, cat, switch):

            if cat is None:
                self.path = os.path.join(self.root, switch)
            else:
                if switch is None:
                    self.path = os.path.join(self.root, cat)
                else:
                    self.path = os.path.join(self.root, cat, switch)

            # self.path = os.path.abspath(self.path)

            # self.files = [os.path.join(dp, f)
            #               for dp, dn, fn in os.walk(self.path)
            #               for f in fn
            #               if f.endswith('F90')
            #               ]

            self.files = [os.path.join(self.path, f)
                          for f in os.listdir(self.path)
                          if f.endswith('F90') or f.endswith('h')]
            self.name = switch

    def __init__(self, root_path, dim):

        # dim specifies root directory.
        self.root = os.path.join(root_path, dim)
        self.item.root = self.root

        # essential things
        self.Driver = self.item(None, 'Driver')
        self.IO = self.item(None, 'IO')
        self.Simulation = self.item(None, 'Simulation')
        self.Numerics = self.item(None, 'Numerics')

    def push_item(self, cat, switch):
        setattr(self, cat.rsplit('/')[-1], self.item(cat, switch))

    def make_symlinks(self):
        for attr, value in self.__dict__.items():
            try:
                fl = getattr(value, 'files')
                for f in fl:
                    dest = os.path.basename(f)
                    try:
                        os.symlink(f, dest)
                    except FileExistsError:
                        print("Removing exsisting symlink: ", dest)
                        os.remove(dest)
                        os.symlink(f, dest)
            except AttributeError:
                pass


def FDM2D(args=None):

    src_path = '../src'
    dimension = '2D'

    # init config
    config = slug_config(src_path, dimension)

    # build the slugcode!
    config.push_item('System', 'Hydro')
    config.push_item('Mesh', 'UniformGrid')
    # I guess temporal doesn't need to have a switch
    config.push_item('Numerics/Temporal', None)
    config.push_item('Numerics/Spatial', 'FDM')

    # copy example (vortex) slug.init file
    if not os.path.exists('./slug.init'):
        shutil.copy(os.path.join(config.root, 'slug.init'), './')

    # link makefile
    if not os.path.exists('./Makefile'):
        os.symlink(os.path.join(config.root, 'Makefile'), './Makefile')
        os.symlink(os.path.join(src_path, 'Makefile.header'), './Makefile.header')

    # print(config.Mesh.files)
    # print(config.IO.files)
    # print(config.System.files)
    # print(config.Numerics.files)
    # print(config.Temporal.files)

    config.make_symlinks()


if __name__ == '__main__':

    FDM2D()
