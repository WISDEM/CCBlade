#!/usr/bin/env python
# encoding: utf-8

# setup.py
# only if building in place: ``python setup.py build_ext --inplace``
import os
import re
import platform
import shutil
import setuptools
import subprocess

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


def run_meson_build(staging_dir):
    prefix = os.path.join(os.getcwd(), staging_dir)
    purelibdir = "."

    # check if meson extra args are specified
    meson_args = ""
    if "MESON_ARGS" in os.environ:
        meson_args = os.environ["MESON_ARGS"]
        # A weird add-on on mac github action runners needs to be removed
        if meson_args.find("buildtype") >= 0: meson_args = ""

    if platform.system() == "Windows":
        if not "FC" in os.environ:
            os.environ["FC"] = "gfortran"
        if not "CC" in os.environ:
            os.environ["CC"] = "gcc"

    # configure
    meson_path = shutil.which("meson")
    if meson_path is None:
        raise OSError("The meson command cannot be found on the system")
        
    meson_call = [meson_path, "setup", staging_dir, "--wipe",
                  f"--prefix={prefix}", f"-Dpython.purelibdir={purelibdir}",
                  f"-Dpython.platlibdir={purelibdir}", meson_args]
    meson_call = [m for m in meson_call if m != ""]
    print(meson_call)
    p1 = subprocess.run(meson_call, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    os.makedirs(staging_dir, exist_ok=True)
    setup_log = os.path.join(staging_dir, "setup.log")
    with open(setup_log, "wb") as f:
        f.write(p1.stdout)
    if p1.returncode != 0:
        with open(setup_log, "r") as f:
            print(f.read())
        raise OSError(meson_call, f"The meson setup command failed! Check the log at {setup_log} for more information.")

    # build
    meson_call = [meson_path, "compile", "-vC", staging_dir]
    meson_call = [m for m in meson_call if m != ""]
    print(meson_call)
    p2 = subprocess.run(meson_call, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    compile_log = os.path.join(staging_dir, "compile.log")
    with open(compile_log, "wb") as f:
        f.write(p2.stdout)
    if p2.returncode != 0:
        with open(compile_log, "r") as f:
            print(f.read())
        raise OSError(meson_call, f"The meson compile command failed! Check the log at {compile_log} for more information.")


def copy_shared_libraries():
    build_path = os.path.join(staging_dir, "ccblade")
    for root, _dirs, files in os.walk(build_path):
        for file in files:
            # move ccblade to just under staging_dir
            if file.endswith((".so", ".lib", ".pyd", ".pdb", ".dylib", ".dll", ".mod")):
                if ".so.p" in root or ".pyd.p" in root:  # excludes intermediate object files
                    continue
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(staging_dir, new_path)
                new_path = new_path[match.span()[1] + 1 :]
                print(f"Copying build file {file_path} -> {new_path}")
                shutil.move(file_path, new_path)


if __name__ == "__main__":
    # This is where the meson build system will install to, it is then
    # used as the sources for setuptools
    staging_dir = "meson_build"

    # this keeps the meson build system from running more than once
    if "dist" not in str(os.path.abspath(__file__)):
        cwd = os.getcwd()
        run_meson_build(staging_dir)
        os.chdir(cwd)
        copy_shared_libraries()

    init_file = os.path.join("ccblade", "__init__.py")
    #__version__ = re.findall(
    #    r"""__version__ = ["']+([0-9\.]*)["']+""",
    #    open(init_file).read(),
    #)[0]

    setuptools.setup(cmdclass={'bdist_wheel': bdist_wheel}, distclass=BinaryDistribution)

#os.environ['NPY_DISTUTILS_APPEND_FLAGS'] = '1'
