#!/usr/bin/env python3
import re
import os
import argparse


MOD_REGEX = re.compile(r"^\s*(?P<unit_type>module(?!\s+procedure)|program)\s*(?P<modname>\w*)",
                       re.IGNORECASE)
USE_REGEX = re.compile(r"""^\s*use
(\s*,\s*intrinsic\s*)?(\s*::\s*|\s+)  # Valid separators between "use" and module name
(?P<moduse>\w*)                       # The module name
\s*(, )?\s*(only)?\s*(:)?.*?$         # Stuff that might follow the name
""",
                       re.IGNORECASE | re.VERBOSE)


def get_filelist(root_path, ext):
    """Get a list of files in 'root_path' with 'ext' extension.

    Args:
        root_path (str): A root directory of Fortran files.
        ext (str): File extension for Fortran files.

    Returns:
        list of filenames (str)
    """

    fl = [os.path.join(root_path, f)
          for f in os.listdir(root_path) if f.endswith(ext)]

    return fl


def get_moddep_dict(file_list):
    """Read Fortran files in given list of files 'file_list',
        and find which modules are needed for a file.

    Args:
        file_list (list of str): A list of Fortran files.

    Returns:
        dictionary of FileName: ModuleName pairs.
    """

    mod_dep = {}
    for filename in file_list:

        with open(filename) as f:
            contents = f.read()
            lines = contents.splitlines()

            uses = []
            for num, line in enumerate(lines):

                found = re.match(USE_REGEX, line)
                if found:
                    uses.append(found.group('moduse').strip())

        mod_dep[os.path.basename(filename)] = uses

    return mod_dep


def get_modpath_dict(file_list):
    """Read Fortran files in given list of files 'file_list',
        and find where the modules are defined at.

    Args:
        file_list (list of str): A list of Fortran files.

    Returns:
        dictionary of ModuleName: ModuleFileName pairs.
    """

    mod_path = {}
    for filename in file_list:
        with open(filename) as f:
            contents = f.read()
            lines = contents.splitlines()

            for num, line in enumerate(lines):

                found = re.match(MOD_REGEX, line)
                if found:
                    mod_path[found.group('modname')] = os.path.basename(filename)

    return mod_path


def get_file_dep_dict(mod_dep_dict, mod_path_dict):
    """Transform ModuleName in 'mod_dep_dict' with ModuleFileName
        with given given 'mod_path_dict'.

    Args:
        mod_dep_dict (dictionary): output of 'get_moddep_dict' function
        mod_path_dict (dictionary): output of 'get_modpath_dict' function

    Returns:
        dictionary of FileName: ModuleFileName pairs.
    """

    file_dep = {}
    warn_printed = []
    for filename, modlist in mod_dep_dict.items():
        mod_file_path = []
        for modname in modlist:
            try:
                mod_file_path.append(mod_path_dict[modname])
            except KeyError:
                if modname not in warn_printed:
                    # suppress multiple warning messages
                    print('[MAKEDEP.PY]: ' + modname +
                          ' is not defined in project folder. '
                          'Could be an external library. PASS')
                    warn_printed.append(modname)
                pass
        if mod_file_path:
            # check if the list is empty
            file_dep[filename] = mod_file_path

    return file_dep



def write_make_dep(file_dep_dict, makefilename):
    """Write file dependency information in given filename.

    Args:
        file_dep_dict (dictionary): file dependency information.
                                    output of 'get_file_dep_dict' function.
        makefilename (str): output file name.

    Returns:
        None
    """

    with open(makefilename, "w") as make:
        header = '# This file is generated automatically'
        make.write(header + '\n')

        for filename, modlist in file_dep_dict.items():
            fileobj = os.path.splitext(filename)[0]+".o"

            modstr = ""
            for modfile in modlist:
                modobj = os.path.splitext(modfile)[0]+".o"
                modstr += " "+str(modobj)

            line = str(fileobj)+":"+modstr+"\n"

            make.write(line)


def create_argument_parser():
    """Create the parser for the command line arguments

    """
    # Add command line arguments
    parser = argparse.ArgumentParser(description='Generate Fortran dependencies')
    # root directory
    parser.add_argument('-r', '--root', nargs='?', const=1, default='./',
                        help='Root directory containing Fortran files.')
    # file extension
    parser.add_argument('-e', '--ext', nargs='?', const=1, default='F90',
                        help='Extension for Fortran files.')
    # output file name
    parser.add_argument('-o', '--output', nargs='?', const=1, default='makefile.dep',
                        help='Output filename. ')

    return parser


def main(args=None):
    """Main function

    """

    # parse command line args.
    parser = create_argument_parser()
    args = parser.parse_args()

    # get project file list
    files = get_filelist(args.root, args.ext)

    # determine which modules are needed for each file
    mod_dep = get_moddep_dict(files)
    # determine where the modules are defined at
    mod_path = get_modpath_dict(files)

    # determine which files are needed for each file
    file_dep = get_file_dep_dict(mod_dep, mod_path)

    # write output
    write_make_dep(file_dep, args.output)


if __name__ == "__main__":

    main()
