#VisIt plugin for SlugCode

It enables [VisIt] to read SlugCode's HDF5 output.
Currently, 1D output is not yet supported.

##Install

Assume that [VisIt] is installed under the `PATH`

```sh
xml2cmake SlugCode.xml
mkdir build
cd build
cmake ..
make
```

Then, libraries will be installed under `~/.visit`.

[VisIt](http://http://visitusers.org) 
