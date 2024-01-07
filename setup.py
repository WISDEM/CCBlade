#!/usr/bin/env python
import os
import re
import shutil
import platform
import setuptools
from setuptools.command.build_ext import build_ext

#######
# This forces wheels to be platform specific
from setuptools.dist import Distribution
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        self.root_is_pure = False

class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name"""
    def has_ext_modules(foo):
        return True
#######
this_dir = os.path.abspath(os.path.dirname(__file__))
staging_dir = os.path.join(this_dir, "meson_build")
build_dir = os.path.join(this_dir, "build")

def copy_shared_libraries():
    build_path = os.path.join(build_dir, "ccblade")
    for root, _dirs, files in os.walk(build_path):
        for file in files:
            if file.endswith((".so", ".lib", ".pyd", ".pdb", ".dylib", ".dll")):
                if ".so.p" in root or ".pyd.p" in root:  # excludes intermediate object files
                    continue
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(build_dir, new_path)
                new_path = new_path[match.span()[1] + 1 :]
                print(f"Copying build file {file_path} -> {new_path}")
                shutil.copy(file_path, new_path)

#######
class MesonExtension(setuptools.Extension):

    def __init__(self, name, sourcedir="", **kwa):
        setuptools.Extension.__init__(self, name, sources=[], **kwa)
        self.sourcedir = os.path.abspath(sourcedir)

class MesonBuildExt(build_ext):
    
    def copy_extensions_to_source(self):
        newext = []
        for ext in self.extensions:
            if isinstance(ext, MesonExtension): continue
            newext.append( ext )
        self.extensions = newext
        super().copy_extensions_to_source()
    
    def build_extension(self, ext):
        if not isinstance(ext, MesonExtension):
            super().build_extension(ext)

        else:

            # Ensure that Meson is present and working
            try:
                self.spawn(["meson", "--version"])
            except OSError:
                raise RuntimeError("Cannot find meson executable")
            
            # check if meson extra args are specified
            meson_args = ""
            if "MESON_ARGS" in os.environ:
                meson_args = os.environ["MESON_ARGS"]

            if platform.system() == "Windows":
                if "FC" not in os.environ:
                    os.environ["FC"] = "gfortran"
                if "CC" not in os.environ:
                    os.environ["CC"] = "gcc"

            purelibdir = "."
            configure_call = ["meson", "setup", staging_dir, "--wipe",
                          f"-Dpython.purelibdir={purelibdir}", f"--prefix={build_dir}", 
                          f"-Dpython.platlibdir={purelibdir}"] + meson_args.split()
            configure_call = [m for m in configure_call if m.strip() != ""]
            print(configure_call)

            build_call = ["meson", "compile", "-vC", staging_dir]
            print(build_call)

            install_call = ["meson", "install", "-C", staging_dir]
            print(install_call)

            self.build_temp = build_dir

            self.spawn(configure_call)
            self.spawn(build_call)
            self.spawn(install_call)
            copy_shared_libraries()

            
if __name__ == "__main__":
    setuptools.setup(cmdclass={"bdist_wheel": bdist_wheel, "build_ext": MesonBuildExt},
                     distclass=BinaryDistribution,
                     ext_modules=[ MesonExtension("ccblade", this_dir) ],
                     )
