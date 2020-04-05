# pylint: disable-msg=W0614,W0401,W0611,W0622
# flake8: noqa

dependencies = ("numpy", "matplotlib")
missing_dependencies = []

for dependency in dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(dependency)

if missing_dependencies:
    raise ImportError(
        "Missing required dependencies {0}".format(missing_dependencies))
del dependencies, dependency, missing_dependencies


from slugger.core.data import *
from slugger.core.pars import *
