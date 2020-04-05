def pars2dict(init_file):
    """ It reads parameters `init_file` and returns a dictionary """

    init_dict = {}

    with open(init_file) as f:
        lines = f.readlines()

        for line in lines:

            # if the line start with #, it is a comment
            if line[0] == "#":
                continue

            # check the line contains a comment
            comment = line.find("#")

            # strip the comment
            if comment != -1:
                line = line[0:comment]

            # skip the empty line
            if line == "\n":
                continue

            # stored :: list
            stored = line.rstrip("\n").strip().split()

            # update stored value as a dictionary
            init_dict[stored[0]] = stored[1]

    return init_dict


def dict2file(init_dict):
    """ Take parameters as a dictionary and write it as a file `slug.init`. """
    f = open("slug.init", "w")
    for key, value in init_dict.items():
        line = str(key) + " " + str(value) + "\n"
        f.write(line)
