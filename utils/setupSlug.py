#!/usr/bin/env python3
import os
import shutil
import argparse


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


def copy_inits(dim):

    # copy init files
    init_dir = os.path.join('../inits', dim)
    init_files = [f for f in os.listdir(init_dir)
                  if f.endswith('.init')]
    for init in init_files:
        init_orig = os.path.join(init_dir, init)
        init_base = os.path.join('.', init)
        shutil.copy(init_orig, init_base)


def FDM1D(args=None):

    src_path = '../src'
    dimension = '1D'

    # init config
    config = slug_config(src_path, dimension)

    # build the slugcode!
    config.push_item('System', 'Hydro')
    config.push_item('Mesh', 'UniformGrid')
    # I guess temporal doesn't need to have a switch
    config.push_item('Numerics/Temporal', None)
    config.push_item('Numerics/Spatial', 'FDM')

    # copy example slug.init file
    if not os.path.exists('./slug.init'):
        shutil.copy(os.path.join(config.root, 'slug.init'), './')

    # link makefile
    if not os.path.exists('./Makefile'):
        os.symlink(os.path.join(config.root, 'Makefile'), './Makefile')
        os.symlink(os.path.join(src_path, 'Makefile.header'), './Makefile.header')

    # copy init files
    copy_inits(dimension)

    # symlinks!
    config.make_symlinks()


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
    config.push_item('Numerics/Positivity', 'MPP')

    # copy example (vortex) slug.init file
    if not os.path.exists('./slug.init'):
        shutil.copy(os.path.join(config.root, 'slug.init'), './')

    # link makefile
    if not os.path.exists('./Makefile'):
        os.symlink(os.path.join(config.root, 'Makefile'), './Makefile')
        os.symlink(os.path.join(src_path, 'Makefile.header'), './Makefile.header')

    # copy init files
    copy_inits(dimension)

    # print(config.Mesh.files)
    # print(config.IO.files)
    # print(config.System.files)
    # print(config.Numerics.files)
    # print(config.Temporal.files)

    # symlinks!
    config.make_symlinks()


def FDM3D(args=None):

    src_path = '../src'
    dimension = '3D'

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

    # copy init files
    copy_inits(dimension)

    # symlinks!
    config.make_symlinks()


def write_version_header():

    commit = os.popen('git rev-parse HEAD').read()
    line = "#define VERSION " + f'"{commit.strip()}"'

    with open('git_version.h', 'w') as f:
        f.write('! generated by setupSlug.py\n')
        f.write(line)


def create_argument_parser():

    # Add command line arguments
    parser = argparse.ArgumentParser(
            description='You never fail until you stop trying.')
    # FDM1D
    parser.add_argument('--fdm1d', action='store_true',
                        help='FDM/1D sovling 1D stencils')
    # FDM2D
    parser.add_argument('--fdm2d', action='store_true',
                        help='FDM/2D sovling 1D stencils')
    # FDM3D
    parser.add_argument('--fdm3d', action='store_true',
                        help='FDM/3D sovling 1D stencils')

    return parser


if __name__ == '__main__':

    parser = create_argument_parser()
    args = parser.parse_args()

    if args.fdm1d:
        FDM1D()
    elif args.fdm2d:
        FDM2D()
    elif args.fdm3d:
        FDM3D()
    else:
        print('Please feed me!')
        parser.print_help()

    write_version_header()
